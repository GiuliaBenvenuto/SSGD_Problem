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
#include <unordered_map>
#include <string>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <Eigen/SparseCholesky>

// Compute SSGD
#include "solving_ssgd.h"

// Matlab
#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"

using namespace std;
using namespace cinolib;
using namespace gcHeatWrapper;
using namespace matlab::engine;
namespace fs = std::filesystem;
using namespace std::chrono;

struct State {

    DrawableTrimesh<> m;
    uint nverts; 

    // Solver
    VTPSolver           vtp_solver;
    // TrettnerSolver      trettner_solver; 
    // FastMarchingSolver  fast_mar_solver;
    // HeatSolver          heat_solver;
    // GeotangleSolver     geotangle_solver;
    // EdgeSolver          edge_solver;
    // ExtendedSolver      extended_solver;
    // LanthierSolver      lanthier_solver;

    // Timer
    double vtp_load,        vtp_preprocess,         vtp_query;
    double trettner_load,   trettner_preprocess,    trettner_query;
    double fast_mar_load,   fast_mar_preprocess,    fast_mar_query;
    double heat_load,       heat_preprocess,        heat_query;
    double geotangle_load,  geotangle_preprocess,   geotangle_query;
    double edge_load,       edge_preprocess,        edge_query;
    double extended_load,   extended_preprocess,    extended_query;
    double lanthier_load,   lanthier_preprocess,    lanthier_query;

    vector<double> blub_ground_truth, bob_ground_truth, spot_ground_truth;

    string mesh_path_1;
    string mesh_path_2;

    vector<double> res;

    State() {
        m = DrawableTrimesh<>();
        nverts = 0;

        blub_ground_truth = vector<double>();
        bob_ground_truth = vector<double>();
        spot_ground_truth = vector<double>();

        mesh_path_1 = "";
        mesh_path_2 = "";

        res = vector<double>();

        vtp_solver = VTPSolver();
    }
};


// Function to read the ground truth for a specific vertex
bool read_ground_truth(const string& csv_gt, int vertex, vector<double>& ground_truth) {
    ifstream gt_csvFile(csv_gt);
    if (!gt_csvFile.is_open()) {
        cerr << "Failed to open ground truth file: " << csv_gt << endl;
        return false;
    }

    string line;
    getline(gt_csvFile, line); 
    stringstream header(line);
    string headerValue;

    int columnIndex = -1;
    int index = 0;

    // Find the correct column for the vertex
    while (getline(header, headerValue, ',')) {
        if (headerValue == "vertex_" + to_string(vertex)) {
            columnIndex = index;
            break;
        }
        index++;
    }

    if (columnIndex == -1) {
        cerr << "Vertex " << vertex << " not found in header." << endl;
        return false;
    } else {
        cout << "Vertex " << vertex << " found at column " << columnIndex << endl;
    }

    // Read the specified column from each line
    while (getline(gt_csvFile, line)) {
        stringstream lineStream(line);
        string cell;
        int currentColumn = 0;

        while (getline(lineStream, cell, ',')) {
            if (currentColumn == columnIndex) {
                ground_truth.push_back(stod(cell));
                break;
            }
            currentColumn++;
        }
    }

    if (ground_truth.empty()) {
        cerr << "Ground truth is empty." << endl;
        return false;
    } else {
        cout << ground_truth.size() << " ground truth values read." << endl;
    }

    gt_csvFile.close();
    return true;
}


// Load meshes
void load_mesh(const string &filename, State &gs) {
    // clear the mesh
    gs.m.clear();
    // clear number of vertices
    gs.nverts = 0;

    gs.m = DrawableTrimesh<>(filename.c_str());
    gs.nverts = gs.m.num_verts();
}

void init(GeodesicMethod &method, State &gs, const string &name) {
    cout << "LOADING " << name << " ..." << endl;
    method.load(&gs.m);
    cout << "PREPROCESSING " << name << " ..." << endl;
    method.preprocess();
}

void init_methods(State &gs) {
    init(gs.vtp_solver, gs, "VTP");
}

void run_ssgd_method(State &gs, int vertex) {
    gs.res.clear();
    cout << "Running SSGE method..." << endl;
    gs.vtp_solver.query(vertex, gs.res);
    cout << "SSGE method completed." << endl;
    gs.res.resize(gs.nverts);

    for (int i = 0; i < 20; i++) {
        cout << "Index " << i << ": " << gs.res[i] << endl;
    }
}


// PRINT VTP bob_tri_final.obj - bob_tri_subdiv_3_final.obj for vertex 482
// Assume other necessary headers and namespace declarations are here
int main(int argc, char **argv) {
    State gs;

    // Valid vertices for the meshes
    vector<int> vv_bob = {482};

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <folder_path>" << endl;
        return 1;
    }
    // folderPath = ../pymeshlab/Esperimento_1/data/bob
    string folderPath = argv[1];

    gs.mesh_path_1 = folderPath + "/bob_tri_final.obj";
    gs.mesh_path_2 = folderPath + "/bob_tri_subdiv_3_final.obj";

    for (int vertex : vv_bob) {
        cout << "MESH: " << gs.mesh_path_1 << endl;
        load_mesh(gs.mesh_path_1, gs);
        init_methods(gs);
        run_ssgd_method(gs, vertex);

        cout << "MESH: " << gs.mesh_path_2 << endl;
        load_mesh(gs.mesh_path_2, gs);
        init_methods(gs);
        run_ssgd_method(gs, vertex);

    }

    return 0;
}


// BACKUP
// int main(int argc, char **argv) {
//     State gs;

//     // Valid vertices for the meshes
//     vector<int> vv_bob = {482, 1710, 2005, 3782, 4757};

//     if (argc < 2) {
//         cerr << "Usage: " << argv[0] << " <folder_path>" << endl;
//         return 1;
//     }
//     string folderPath = argv[1];

//     for (int vertex : vv_bob) {
//         cout << "------- Processing vertex: " << vertex << " --------" << endl;
//         gs.bob_ground_truth.clear();

//         // read a csv file
//         string csv_gt = "../pymeshlab/Esperimento_1/data/gt/bob_gt_prova.csv";
//         cout << "Reading ground truth from: " << csv_gt << endl;

//         if (!read_ground_truth(csv_gt, vertex, gs.bob_ground_truth)) {
//             return 1;
//         } else {
//             cout << "Ground truth read successfully." << endl;
//             // for (int i = 0; i < 10; i++) {
//             //     cout << "Index " << i << ": " << gs.bob_ground_truth[i] << endl;
//             // }
//         }
//     }

//     return 0;
// }