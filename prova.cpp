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
#include <cmath>

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
    GeotangleSolver     geotangle_solver;
    // EdgeSolver          edge_solver;
    ExtendedSolver      extended_solver;
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

    // string mesh_path_1, mesh_path_2;
    string mesh_path;
    string mesh_name;

    vector<double> res;

    State() {
        m = DrawableTrimesh<>();
        nverts = 0;

        blub_ground_truth = vector<double>();
        bob_ground_truth = vector<double>();
        spot_ground_truth = vector<double>();

        // mesh_path_1 = "";
        // mesh_path_2 = "";
        mesh_path = "";
        mesh_name = "";

        res = vector<double>();

        // vtp_solver = VTPSolver();
        // geotangle_solver = GeotangleSolver();
    }
};


// Function to read the ground truth for a specific vertex
// HO VERIFICATO CHE FUNZIONA CORRETTAMENTE
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
    cout << "load mesh function" << endl;
    // Reset the mesh
    gs.m.clear();
    gs.nverts = 0;

    gs.m = DrawableTrimesh<>(filename.c_str());
    gs.nverts = gs.m.num_verts();
    cout << "Mesh loaded: " << filename << endl;
    cout << "Number of vertices: " << gs.nverts << endl;
}

void init(GeodesicMethod &method, State &gs, const string &name) {
    cout << "LOADING " << name << " ..." << endl;
    method.load(&gs.m);
    cout << "PREPROCESSING " << name << " ..." << endl;
    method.preprocess();
}

void init_methods(State &gs) {
    // init(gs.vtp_solver, gs, "VTP");

    // gs.geotangle_solver = GeotangleSolver();
    // init(gs.geotangle_solver, gs, "Geotangle");

    init(gs.extended_solver, gs, "Extended");

    
}


double smape_error(const vector<double> &res, const vector<double> &gt, State &gs) {
    // gt: ground truth
    // res: result obtained from the method
    cout << "Res size: " << res.size() << endl;
    cout << "GT size: " << gt.size() << endl;

    if (res.size() != gs.nverts) {
        cerr << "Error: res and nverts have different sizes." << endl;
        return -1;
    }

    // if(res.size() != gt.size()) {
    //     cerr << "Error: res and gt have different sizes." << endl;
    //     return -1;
    // }
    double smape = 0.0;
    int count = 0;

    for (int i = 0; i < res.size(); i++) {
        double num = abs(res[i] - gt[i]);
        // double den = (abs(res[i]) + abs(gt[i])) / 2.0;
        double den = (abs(res[i]) + abs(gt[i]));

        if (den != 0) {
            smape += num / den;
            count++;
        }
    }

    if (count == 0) {
        cerr << "Error: No valid data points to calculate SMAPE." << endl;
        return -1; // Or handle this case differently if preferred
    }
    cout << "Smape: " << smape << endl;
    cout << "Count: " << count << endl;

    return (smape / count) * 100.0;
}


// void run_ssgd_method(State &gs, int vertex, ofstream &smape_csv, ofstream &dist_csv) {
void run_ssgd_method(State &gs, int vertex) {
    // gs.res.clear();
    // cout << "Running SSGE method..." << endl;
    // gs.vtp_solver.query(vertex, gs.res);
    // cout << "SSGE method completed." << endl;
    // gs.res.resize(gs.nverts);

    gs.res.clear();
    cout << "Running Extended method..." << endl;
    gs.extended_solver.query(vertex, gs.res);
    cout << "Extended method completed." << endl;
    gs.res.resize(gs.nverts);


    // gs.res.clear();
    // cout << "Running Geotangle method..." << endl;
    // gs.geotangle_solver.query(vertex, gs.res);
    // cout << "Geotangle method completed." << endl;
    // gs.res.resize(gs.nverts);

    // // call the smape error function here
    // double smape = smape_error(gs.res, gs.bob_ground_truth, gs);
    // cout << "SMAPE: " << smape << endl;

    // // Write the results to the csv file
    // smape_csv << gs.mesh_name << "," << gs.nverts << "," << vertex << "," << smape << endl;

    // for (int i = 0; i < gs.res.size(); i++) {
    //     dist_csv << i << "," << gs.res[i] << endl;
    // }
}


// BACKUP
int main(int argc, char **argv) {
    State gs;

    // Valid vertices for the meshes
    // vector<int> vv_bob = {482, 1710, 2005, 3782, 4757};
    // vector<int> vv_bob = {482, 1710, 2005, 3782, 4757};
    vector<int> vv_bob = {482};

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <folder_path>" << endl;
        return 1;
    }
    string folderPath = argv[1];

    // OUTPUT CSV
    ofstream smape_csv("../pymeshlab/Esperimento_1/data/prova/SMAPE_bob.csv");
    if (!smape_csv.is_open()) {
        cerr << "Failed to open SMAPE file." << endl;
        return 1;
    }
    smape_csv << "Mesh,Vertices,Vertex,SMAPE_Geotangle" << endl;

    ofstream dist_csv("../pymeshlab/Esperimento_1/data/prova/dist_BOB_ORIGINAL.csv");
    if (!dist_csv.is_open()) {
        cerr << "Failed to open SMAPE file." << endl;
        return 1;
    }
    dist_csv << "Index,Dist_482" << endl;



    // For each vertex
    for (int vertex : vv_bob) {
        cout << endl << endl << "------- Processing vertex: " << vertex << " --------" << endl;
        gs.bob_ground_truth.clear();

        // read a csv file
        string csv_gt = "../pymeshlab/Esperimento_1/data/gt/bob_gt_distances.csv";
        cout << "Reading ground truth from: " << csv_gt << endl;

        if (!read_ground_truth(csv_gt, vertex, gs.bob_ground_truth)) {
            return 1;
        } else {
            cout << "Ground truth read successfully." << endl << endl;
        }

        // For element in the folder with extension ".obj"
        for (const auto &entry : fs::directory_iterator(folderPath)) {
            if (entry.path().extension() == ".obj") {
                cout << endl << "MESH: " << entry.path() << endl;
                gs.mesh_path = entry.path().string();
                gs.mesh_name = entry.path().filename().string();

                load_mesh(gs.mesh_path, gs);
                init_methods(gs);
                // run_ssgd_method(gs, vertex, smape_csv, dist_csv);

            }
        }
    }

    return 0;
}


// // ----- Compute again VTP for subdiv_5 for vertex 482 -----
// int main(int argc, char **argv) {
//     State gs;

//     // Valid vertices for the meshes
//     // vector<int> vv_bob = {482};
//     vector<int> vv_bunny = {100};

//     if (argc < 2) {
//         cerr << "Usage: " << argv[0] << " <folder_path>" << endl;
//         return 1;
//     }
//     // folderPath = ../pymeshlab/Esperimento_1/data/bob
//     string folderPath = argv[1];

//     // gs.mesh_path = folderPath + "/bob_tri_subdiv_5_final.obj";
//     // gs.mesh_path = folderPath + "/bunny_500_subdiv_6_final.obj";
//     gs.mesh_path = folderPath + "/bunny_500_final.obj";

//     for (int vertex : vv_bunny) {
//         cout << "MESH: " << gs.mesh_path << endl;
//         load_mesh(gs.mesh_path, gs);
//         init_methods(gs);
//         run_ssgd_method(gs, vertex);

//         // take gs.res and write it to a csv file in the same folder
//         // string csv_res = "../pymeshlab/Esperimento_1/data/prova/bob_tri_subdiv_5_final_res.csv";
//         string csv_res = "../pymeshlab/Esperimento_1/data/prova/bunny_500_ext_res.csv";
//         ofstream res_csvFile(csv_res);
//         if (!res_csvFile.is_open()) {
//             cerr << "Failed to open result file: " << csv_res << endl;
//             return 1;
//         }
//         res_csvFile << "vertex_" << vertex << endl;
//         for (int i = 0; i < gs.res.size(); i++) {
//             res_csvFile << gs.res[i] << endl;
//         }
//         res_csvFile.close();
//     }
//     return 0;
// }



// ----- PRINT VTP bob_tri_final.obj - bob_tri_subdiv_3_final.obj for vertex 482 -----
// Assume other necessary headers and namespace declarations are here

// int main(int argc, char **argv) {
//     State gs;

//     // Valid vertices for the meshes
//     vector<int> vv_bob = {482};

//     if (argc < 2) {
//         cerr << "Usage: " << argv[0] << " <folder_path>" << endl;
//         return 1;
//     }
//     // folderPath = ../pymeshlab/Esperimento_1/data/bob
//     string folderPath = argv[1];

//     gs.mesh_path_1 = folderPath + "/bob_tri_final.obj";
//     gs.mesh_path_2 = folderPath + "/bob_tri_subdiv_3_final.obj";

//     for (int vertex : vv_bob) {
//         cout << "MESH: " << gs.mesh_path_1 << endl;
//         load_mesh(gs.mesh_path_1, gs);
//         init_methods(gs);
//         run_ssgd_method(gs, vertex);

//         cout << "MESH: " << gs.mesh_path_2 << endl;
//         load_mesh(gs.mesh_path_2, gs);
//         init_methods(gs);
//         run_ssgd_method(gs, vertex);

//     }

//     return 0;
// }





