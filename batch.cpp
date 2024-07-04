/*
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
#include "SSGD_methods/heat/gc_wrapper.h"

// Matlab
#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"


using namespace std;
using namespace cinolib;
using namespace gcHeatWrapper;
using namespace matlab::engine;

struct State {
    DrawableTrimesh<> m;    // the input mesh
    uint nverts;            // number of vertices
    vector<uint> tris;      // triangle indices
    vector<double> coords;  // vertex coordinates
    vector<double> res;     // result

    // Cache for Heat method
    GeodesicsCache prefactored_matrices;

    // K for extended
    int k;
    int prev_k;

    // Time for heat
    double heat_time;
    double heat_time_prev;

    // Number of Steiner points
    int n_steiner;
    int prev_n_steiner;

    // Sources for SSGD
    std::vector<uint> sources_heat;
    vector<int> sources;

    // Trettner
    string mesh_path;
    string mesh_name;
    HalfEdge mesh;

    // Solver
    VTPSolver           vtp_solver;
    TrettnerSolver      trettner_solver; 
    FastMarchingSolver  fast_mar_solver;
    HeatSolver          heat_solver;
    GeotangleSolver     geotangle_solver;
    EdgeSolver          edge_solver;
    ExtendedSolver      extended_solver;
    LanthierSolver      lanthier_solver;

    // Timer
    std::chrono::steady_clock::time_point tic;
    std::chrono::steady_clock::time_point toc;
    double vtp_load,        vtp_preprocess,         vtp_query;
    double trettner_load,   trettner_preprocess,    trettner_query;
    double fast_mar_load,   fast_mar_preprocess,    fast_mar_query;
    double heat_load,       heat_preprocess,        heat_query;
    double geotangle_load,  geotangle_preprocess,   geotangle_query;
    double edge_load,       edge_preprocess,        edge_query;
    double extended_load,   extended_preprocess,    extended_query;
    double lanthier_load,   lanthier_preprocess,    lanthier_query;

    State() {
        // Timer
        vtp_load,       vtp_preprocess,       vtp_query = 0.0;
        trettner_load,  trettner_preprocess,  trettner_query = 0.0;
        fast_mar_load,  fast_mar_preprocess,  fast_mar_query = 0.0;
        heat_load,      heat_preprocess,      heat_query = 0.0;
        geotangle_load, geotangle_preprocess, geotangle_query = 0.0;
        edge_load,      edge_preprocess,       edge_query = 0.0;
        extended_load,  extended_preprocess,  extended_query = 0.0;
        lanthier_load,  lanthier_preprocess,  lanthier_query = 0.0;

        res = vector<double>();

        mesh_path = "";
        mesh_name = "";
        mesh = HalfEdge();

        k = 3;
        prev_k = 3;

        n_steiner = 3;
        prev_n_steiner = 3;

        heat_time = 1.0;
        heat_time_prev = 1.0;
    }
};

// ------ Helper functions ------
vector<double> extract_coords(const DrawableTrimesh<> &mesh) {
  vector<double> coords;
  auto verts = mesh.vector_verts();
  for (const auto& vert : verts) {
      coords.push_back(vert.x());
      coords.push_back(vert.y());
      coords.push_back(vert.z());
  }
  return coords;
}

vector<uint> extract_tris(const DrawableTrimesh<> &mesh) {
  vector<uint> tris;
  auto polys = mesh.vector_polys();
  for (const auto& poly : polys) {
      for (auto vid : poly) {
          tris.push_back(vid);
      }
  }
  return tris;
}

void load_mesh(const string &filename, State &gs) {
    gs.m = DrawableTrimesh<>(filename.c_str());
    gs.nverts = gs.m.num_verts();
    gs.coords = extract_coords(gs.m);
    gs.tris = extract_tris(gs.m);
    gs.trettner_solver.mesh_path = filename;
}


void run_ssgd_method(State &state, int sourceVertexIndex) {
    vector<double> distances;
    ScalarField field;
    string resultsFilename = "results.csv";  // Consider adding a timestamp or an identifier to make filenames unique if needed

    ofstream csvFile(resultsFilename, ios::trunc); // Truncate mode to overwrite data for each run
    csvFile << "Mesh Name,Total Vertices,Source Vertex Index,Method,Step,Time (ms),Target Vertex Index,Distance\n";

    auto log_time = [&](auto &solver, const string &method, const string &step) {
        auto tic = chrono::high_resolution_clock::now();

        if (step == "Load") solver.load(state.coords, state.tris);
        else if (step == "Preprocess") solver.preprocess();
        else if (step == "Query") solver.query(sourceVertexIndex, distances, field);

        auto toc = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(toc - tic).count();
        double time_ms = duration / 1000.0;
        cout << step << " time: " << time_ms << " ms" << endl;

        csvFile << state.mesh_name << "," << state.nverts << "," << sourceVertexIndex << ","
                << method << "," << step << "," << time_ms << ",,\n";
    };


    // Helper lambda to log distances
    auto log_distances = [&](const string &method) {
        for (size_t i = 0; i < distances.size(); ++i) {
            csvFile << state.mesh_name << "," << state.nverts << "," << sourceVertexIndex << ","
                    << method << ",,," << i << "," << distances[i] << "\n";
        }
    };


    // VTP Solver
    cout << endl << "----- VTP -----" << endl;
    log_time(state.vtp_solver, "VTP", "Load");
    log_time(state.vtp_solver, "VTP", "Preprocess");
    log_time(state.vtp_solver, "VTP", "Query");
    log_distances("VTP");

    // Trettner Solver
    cout << endl << "----- Trettner -----" << endl;
    log_time(state.trettner_solver, "Trettner", "Load");
    log_time(state.trettner_solver, "Trettner", "Preprocess");
    log_time(state.trettner_solver, "Trettner", "Query");
    log_distances("Trettner");

    // Fast Marching Solver
    cout << endl << "----- Fast Marching -----" << endl;
    log_time(state.fast_mar_solver, "Fast Marching", "Load");
    log_time(state.fast_mar_solver, "Fast Marching", "Preprocess");
    log_time(state.fast_mar_solver, "Fast Marching", "Query");
    log_distances("Fast Marching");

    // Heat Solver
    cout << endl << "----- Heat -----" << endl;
    log_time(state.heat_solver, "Heat", "Load");
    log_time(state.heat_solver, "Heat", "Preprocess");
    log_time(state.heat_solver, "Heat", "Query");
    log_distances("Heat");

    // Geotangle Solver
    cout << endl << "----- Geotangle -----" << endl;
    log_time(state.geotangle_solver, "Geotangle", "Load");
    log_time(state.geotangle_solver, "Geotangle", "Preprocess");
    log_time(state.geotangle_solver, "Geotangle", "Query");
    log_distances("Geotangle");

    // Edge Solver
    cout << endl << "----- Edge -----" << endl;
    log_time(state.edge_solver, "Edge", "Load");
    log_time(state.edge_solver, "Edge", "Preprocess");
    log_time(state.edge_solver, "Edge", "Query");
    log_distances("Edge");

    // Extended Solver
    cout << endl << "----- Extended -----" << endl;
    log_time(state.extended_solver, "Extended", "Load");
    log_time(state.extended_solver, "Extended", "Preprocess");
    log_time(state.extended_solver, "Extended", "Query");
    log_distances("Extended");

    // Lanthier Solver
    cout << endl << "----- Lanthier -----" << endl;
    log_time(state.lanthier_solver, "Lanthier", "Load");
    log_time(state.lanthier_solver, "Lanthier", "Preprocess");
    log_time(state.lanthier_solver, "Lanthier", "Query");
    log_distances("Lanthier");


    csvFile.close();
}

int main(int argc, char **argv) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <mesh_path> <source_vertex_index>" << endl;
        return 1;
    }

    string meshPath = argv[1];
    int sourceVertexIndex = stoi(argv[2]);

    State gs;
    gs.mesh_path = meshPath;
    gs.mesh_name = meshPath.substr(meshPath.find_last_of("/\\") + 1);
    load_mesh(meshPath, gs);
    run_ssgd_method(gs, sourceVertexIndex);

    return 0;
}
*/


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
#include "SSGD_methods/heat/gc_wrapper.h"

// Matlab
#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"

using namespace std;
using namespace cinolib;
using namespace gcHeatWrapper;
using namespace matlab::engine;
namespace fs = std::filesystem;

struct State {
    DrawableTrimesh<> m;    // the input mesh
    uint nverts;            // number of vertices
    vector<uint> tris;      // triangle indices
    vector<double> coords;  // vertex coordinates
    vector<double> res;     // result

    // Cache for Heat method
    GeodesicsCache prefactored_matrices;

    // K for extended
    int k;
    int prev_k;

    // Time for heat
    double heat_time;
    double heat_time_prev;

    // Number of Steiner points
    int n_steiner;
    int prev_n_steiner;

    // Sources for SSGD
    std::vector<uint> sources_heat;
    vector<int> sources;

    // Trettner
    string mesh_path;
    string mesh_name;
    HalfEdge mesh;

    // Solver
    VTPSolver           vtp_solver;
    TrettnerSolver      trettner_solver; 
    FastMarchingSolver  fast_mar_solver;
    HeatSolver          heat_solver;
    GeotangleSolver     geotangle_solver;
    EdgeSolver          edge_solver;
    ExtendedSolver      extended_solver;
    LanthierSolver      lanthier_solver;

    // Timer
    std::chrono::steady_clock::time_point tic;
    std::chrono::steady_clock::time_point toc;
    double vtp_load,        vtp_preprocess,         vtp_query;
    double trettner_load,   trettner_preprocess,    trettner_query;
    double fast_mar_load,   fast_mar_preprocess,    fast_mar_query;
    double heat_load,       heat_preprocess,        heat_query;
    double geotangle_load,  geotangle_preprocess,   geotangle_query;
    double edge_load,       edge_preprocess,        edge_query;
    double extended_load,   extended_preprocess,    extended_query;
    double lanthier_load,   lanthier_preprocess,    lanthier_query;

    State() {
        // Timer
        vtp_load = vtp_preprocess = vtp_query = 0.0;
        trettner_load = trettner_preprocess = trettner_query = 0.0;
        fast_mar_load = fast_mar_preprocess = fast_mar_query = 0.0;
        heat_load = heat_preprocess = heat_query = 0.0;
        geotangle_load = geotangle_preprocess = geotangle_query = 0.0;
        edge_load = edge_preprocess = edge_query = 0.0;
        extended_load = extended_preprocess = extended_query = 0.0;
        lanthier_load = lanthier_preprocess = lanthier_query = 0.0;

        res = vector<double>();

        mesh_path = "";
        mesh_name = "";
        mesh = HalfEdge();

        k = 3;
        prev_k = 3;

        n_steiner = 3;
        prev_n_steiner = 3;

        heat_time = 1.0;
        heat_time_prev = 1.0;
    }
};

// ------ Helper functions ------
vector<double> extract_coords(const DrawableTrimesh<> &mesh) {
    vector<double> coords;
    auto verts = mesh.vector_verts();
    for (const auto& vert : verts) {
        coords.push_back(vert.x());
        coords.push_back(vert.y());
        coords.push_back(vert.z());
    }
    return coords;
}

vector<uint> extract_tris(const DrawableTrimesh<> &mesh) {
    vector<uint> tris;
    auto polys = mesh.vector_polys();
    for (const auto& poly : polys) {
        for (auto vid : poly) {
            tris.push_back(vid);
        }
    }
    return tris;
}

void load_mesh(const string &filename, State &gs) {
    gs.m = DrawableTrimesh<>(filename.c_str());
    gs.nverts = gs.m.num_verts();
    gs.coords = extract_coords(gs.m);
    gs.tris = extract_tris(gs.m);
    gs.trettner_solver.mesh_path = filename;
}

void run_ssgd_method(State &state, int sourceVertexIndex, ofstream &csvFile) {
    vector<double> distances;
    ScalarField field;

    auto log_time = [&](auto &solver, const string &method, const string &step) {
        auto tic = chrono::high_resolution_clock::now();

        if (step == "Load") solver.load(state.coords, state.tris);
        else if (step == "Preprocess") solver.preprocess();
        else if (step == "Query") solver.query(sourceVertexIndex, distances, field);

        auto toc = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(toc - tic).count();
        double time_ms = duration / 1000.0;
        cout << step << " time: " << time_ms << " ms" << endl;

        csvFile << state.mesh_name << "," << state.nverts << "," << sourceVertexIndex << ","
                << method << "," << step << "," << time_ms << ",,\n";
    };

    // Helper lambda to log distances
    auto log_distances = [&](const string &method) {
        for (size_t i = 0; i < distances.size(); ++i) {
            csvFile << state.mesh_name << "," << state.nverts << "," << sourceVertexIndex << ","
                    << method << ",,," << i << "," << distances[i] << "\n";
        }
    };

    // VTP Solver
    cout << endl << "----- VTP -----" << endl;
    log_time(state.vtp_solver, "VTP", "Load");
    log_time(state.vtp_solver, "VTP", "Preprocess");
    log_time(state.vtp_solver, "VTP", "Query");
    log_distances("VTP");

    // Trettner Solver
    cout << endl << "----- Trettner -----" << endl;
    log_time(state.trettner_solver, "Trettner", "Load");
    log_time(state.trettner_solver, "Trettner", "Preprocess");
    log_time(state.trettner_solver, "Trettner", "Query");
    log_distances("Trettner");

    // Fast Marching Solver
    cout << endl << "----- Fast Marching -----" << endl;
    log_time(state.fast_mar_solver, "Fast Marching", "Load");
    log_time(state.fast_mar_solver, "Fast Marching", "Preprocess");
    log_time(state.fast_mar_solver, "Fast Marching", "Query");
    log_distances("Fast Marching");

    // Heat Solver
    cout << endl << "----- Heat -----" << endl;
    log_time(state.heat_solver, "Heat", "Load");
    log_time(state.heat_solver, "Heat", "Preprocess");
    log_time(state.heat_solver, "Heat", "Query");
    log_distances("Heat");

    // Geotangle Solver
    cout << endl << "----- Geotangle -----" << endl;
    log_time(state.geotangle_solver, "Geotangle", "Load");
    log_time(state.geotangle_solver, "Geotangle", "Preprocess");
    log_time(state.geotangle_solver, "Geotangle", "Query");
    log_distances("Geotangle");

    // Edge Solver
    cout << endl << "----- Edge -----" << endl;
    log_time(state.edge_solver, "Edge", "Load");
    log_time(state.edge_solver, "Edge", "Preprocess");
    log_time(state.edge_solver, "Edge", "Query");
    log_distances("Edge");

    // Lanthier Solver
    cout << endl << "----- Lanthier -----" << endl;
    log_time(state.lanthier_solver, "Lanthier", "Load");
    log_time(state.lanthier_solver, "Lanthier", "Preprocess");
    log_time(state.lanthier_solver, "Lanthier", "Query");
    log_distances("Lanthier");

    // Extended Solver
    cout << endl << "----- Extended -----" << endl;
    log_time(state.extended_solver, "Extended", "Load");
    log_time(state.extended_solver, "Extended", "Preprocess");
    log_time(state.extended_solver, "Extended", "Query");
    log_distances("Extended");
}

int main(int argc, char **argv) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <folder_path> <source_vertex_index>" << endl;
        return 1;
    }

    string folderPath = argv[1];
    int sourceVertexIndex = stoi(argv[2]);

    // Create a new directory for CSV files
    fs::path currentPath = fs::current_path();
    fs::path csvDirPath = currentPath / "csv_idx_1";

    if (!fs::exists(csvDirPath)) {
        fs::create_directory(csvDirPath);
        cout << "Created directory: " << csvDirPath << endl;
    } else {
        cout << "Directory already exists: " << csvDirPath << endl;
    }

    for (const auto &entry : fs::directory_iterator(folderPath)) {
        if (entry.path().extension() == ".obj") {
            string meshPath = entry.path().string();
            State gs;
            gs.mesh_path = meshPath;
            gs.mesh_name = entry.path().filename().string();
            load_mesh(meshPath, gs);
            
            // remove .obj from gs.mesh_name
            gs.mesh_name = gs.mesh_name.substr(0, gs.mesh_name.size() - 4);
            string resultsFilename = (csvDirPath / (gs.mesh_name + "_idx_1.csv")).string();
            ofstream csvFile(resultsFilename, ios::trunc); // Truncate mode to overwrite data for each run
            csvFile << "Mesh Name,Total Vertices,Source Vertex Index,Method,Step,Time (ms),Target Vertex Index,Distance\n";

            run_ssgd_method(gs, sourceVertexIndex, csvFile);
            
            csvFile.close();
        }
    }

    return 0;
}
