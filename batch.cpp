#include <Eigen/SparseCholesky>
#include <chrono>
#include <cinolib/drawable_segment_soup.h>
#include <cinolib/drawable_sphere.h>
#include <cinolib/drawable_vector_field.h>
#include <cinolib/gl/file_dialog_open.h>
#include <cinolib/gl/file_dialog_save.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gradient.h>
#include <cinolib/io/write_OBJ.h>
#include <cinolib/scalar_field.h>
#include <cinolib/vector_serialization.h>
#include <fstream>
#include <imgui.h>
#include <thread>
#include <future>
#include <filesystem>


// Compute SSGD
#include "solving_ssgd.h"

// Matlab
#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"

using namespace std;
using namespace cinolib;
using namespace gcHeatWrapper;
using namespace matlab::engine;
using namespace std::chrono;
namespace fs = std::filesystem;


//--------- CALCOLO GROUND TRUTH per ogni vertice su ogni mesh ---------

struct State {
    DrawableTrimesh<> m;    // the input mesh
    uint nverts;            // number of vertices
    vector<uint> tris;      // triangle indices
    vector<double> coords;  // vertex coordinates
    vector<double> res;     // result
    ScalarField field;

    // Trettner
    string mesh_path;
    string mesh_name;

    // Solver
    VTPSolver           vtp_solver;

    std::chrono::steady_clock::time_point tic;
    std::chrono::steady_clock::time_point toc;

    State() {
        res = vector<double>();

        mesh_path = "";
        mesh_name = "";
    }
};


// Load meshes
void load_mesh(const string &filename, State &gs) {
    gs.m = DrawableTrimesh<>(filename.c_str());
    gs.nverts = gs.m.num_verts();
    gs.coords = extract_coords(gs.m);
    gs.tris = extract_tris(gs.m);
}

void init(GeodesicMethod &m, State &gs, const string &name) {
    cout << "---------- Initializing method: " << name << " ----------" << endl;
    cout << "LOADING..." << endl;
    m.load(&gs.m);

    cout << "PREPROCESSING..." << endl;
    m.preprocess();
    cout << "-------------------------------------------------" << endl;
}

void init_methods(State &gs) {
    init(gs.vtp_solver, gs, "VTP");
}


void run_ssgd_method(State &state, int sourceVertexIndex, vector<vector<double>> &all_distances) {
    vector<double> distances;

    // compute time for VTP
    auto start = chrono::high_resolution_clock::now();
    state.vtp_solver.query(sourceVertexIndex, distances);
    auto end = chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    cout << "QUERY TIME: " << elapsed.count() << " s" << endl;
    cout << "-------------------------------------------------" << endl;

    all_distances.push_back(distances);
}


void write_csv(const string &filename, const vector<int> &vv_blub, const vector<vector<double>> &all_distances) {
    ofstream csvFile(filename);

    if (!csvFile.is_open()) {
        cerr << "Failed to open CSV file for writing." << endl;
        return;
    }

    // Write header
    csvFile << "index";
    for (int vertex : vv_blub) {
        csvFile << ",vertex_" << vertex;
    }
    csvFile << "\n";

    // Write data
    size_t max_size = 0;
    for (const auto &distances : all_distances) {
        max_size = max(max_size, distances.size());
    }

    for (size_t i = 0; i < max_size; ++i) {
        csvFile << i;
        for (const auto &distances : all_distances) {
            if (i < distances.size()) {
                csvFile << "," << distances[i];
            } else {
                csvFile << ",";
            }
        }
        csvFile << "\n";
        csvFile.flush();  // Ensure data is written to file
    }

    csvFile.close();
}


int main(int argc, char **argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <file_path>" << endl;
        return 1;
    }

    State gs;
    string filePath = argv[1];

    vector<int> vv_blub = {663, 3958, 4662, 4715, 6694};

    fs::path csvDirPath = "../pymeshlab/Esperimento_1/data/blub_gt_csv";
    if (!fs::exists(csvDirPath)) {
        fs::create_directory(csvDirPath);
        cout << "Created directory: " << csvDirPath << endl;
    } else {
        cout << "Directory already exists: " << csvDirPath << endl;
    }

    cout << "---------- Loading Mesh ----------" << endl;
    gs.mesh_path = filePath;
    load_mesh(gs.mesh_path, gs);
    cout << "Number of vertices in mesh: " << gs.nverts << endl;
    init_methods(gs);

    vector<vector<double>> all_distances;
    for (int sourceVertexIndex : vv_blub) {
        cout << "Source vertex: " << sourceVertexIndex;
        run_ssgd_method(gs, sourceVertexIndex, all_distances);
    }

    string csvFile = csvDirPath.string() + "/blub_gt_distances.csv";
    write_csv(csvFile, vv_blub, all_distances);

    return 0;
}