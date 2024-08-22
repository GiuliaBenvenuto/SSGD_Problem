// *********************************************************************
// QUESTO E' PER SCRIVERE IL CSV IN CUI HO VTP SULLA MESH PIU GRANDE
// PROGRAMMA PER SCRIVERE IL CSV DELLE GROUND TRUTH PER PIU VERTICI
// *********************************************************************

#include <Eigen/SparseCholesky>
#include <chrono>
#include <cinolib/drawable_segment_soup.h>
#include <cinolib/drawable_sphere.h>
#include <cinolib/drawable_vector_field.h>
#include <cinolib/io/read_OBJ.h>
#include <cinolib/gradient.h>
#include <cinolib/scalar_field.h>
#include <cinolib/vector_serialization.h>
#include <fstream>
#include <thread>
#include <future>

// Compute SSGD
#include "solving_ssgd.h"

// Matlab
#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"


using namespace std;
using namespace cinolib;
using namespace gcHeatWrapper;
namespace fs = std::filesystem;
using namespace matlab::engine;

struct State {
  DrawableTrimesh<> m;
  uint nverts;
  vector<vector<uint>> VV;
  vector<double> coords;
  vector<uint> tris;
  vector<double> res;

  int k;
  int n_stainer;

  VTPSolver vtp_solver;

  vector<int> sources;
};


void Load_mesh(string filename, State &gs) {
    gs.m = DrawableTrimesh<>(filename.c_str());
    gs.nverts = gs.m.num_verts();
    gs.VV.resize(gs.nverts);

    gs.coords = extract_coords(gs.m);
    gs.tris = extract_tris(gs.m);

    for (auto i = 0; i < gs.nverts; i++)
    gs.VV[i] = gs.m.vert_ordered_verts_link(i);

    // QUIIIII
    gs.m.normalize_bbox();
    gs.m.center_bbox();
}


void init(GeodesicMethod &method, State &gs, const string &name) {
    cout << "----- Init method: " << name << " -----" << endl;
    method.load(&gs.m);
    method.preprocess();
}

void init_methods(State &gs) {
    init(gs.vtp_solver, gs, "VTP");
}


// ------- WRITE CSV GROUND TRUTH FILE BACKUP-------
int main(int argc, char **argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <mesh_file.obj>" << endl;
        return 1;
    }

    State gs;
    string mesh_filename = argv[1];

    // Load mesh
    Load_mesh(mesh_filename, gs);

    // Initialize VTP method
    init_methods(gs);

    // Define the array of vertices to compute VTP queries
    vector<int> vv_thai = {48, 50, 56, 49, 2};

    // Define the path to the CSV file
    std::string csv_path = "../pymeshlab/Esperimento_Thai/ground_truth/thai_decimated_gt.csv";
    
    // Create CSV file    
    std::ofstream csv_file(csv_path);
    if (!csv_file.is_open()) {
        cerr << "Failed to create CSV file: " << csv_path << endl;
        return 1;
    }

    // Write CSV header
    csv_file << "index";
    for (const auto& vertex : vv_thai) {
        csv_file << ",vertex_" << vertex;
    }
    csv_file << endl;

    // Vector to store all results
    vector<vector<double>> all_results;

    // Loop over each source vertex and compute VTP distances
    for (const auto& source : vv_thai) {
        cout << "---> Computing VTP distances for vertex " << source << endl;
        // Compute SSGD using VTP method for the current source
        auto tic = std::chrono::steady_clock::now();
        gs.vtp_solver.query(source, gs.res);
        auto toc = std::chrono::steady_clock::now();
        double query_time = chrono::duration_cast<chrono::milliseconds>(toc - tic).count();

        cout << "VTP Query time for vertex " << source << ": " << query_time << " milliseconds" << endl;

        // Store the results in all_results
        if (all_results.empty()) {
            all_results.resize(gs.res.size(), vector<double>(vv_thai.size()));
        }

        for (size_t i = 0; i < gs.res.size(); ++i) {
            all_results[i][&source - &vv_thai[0]] = gs.res[i];
        }
    }

    // Write results to CSV
    for (size_t i = 0; i < all_results.size(); ++i) {
        csv_file << i;
        for (const auto& distance : all_results[i]) {
            csv_file << "," << distance;
        }
        csv_file << endl;
    }

    csv_file.close();
    cout << "Results written to: " << csv_path << endl;

    return 0;
}