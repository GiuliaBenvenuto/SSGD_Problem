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

        vtp_solver = VTPSolver();
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

    gs.geotangle_solver = GeotangleSolver();
    init(gs.geotangle_solver, gs, "Geotangle");
}


double smape_error(const vector<double> &res, const vector<double> &gt, State &gs) {
    // gt: ground truth
    // res: result obtained from the method
    cout << "Res size: " << res.size() << endl;
    cout << "GT size: " << gt.size() << endl;

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

    double smape_error = (smape / count) * 100.0;
    cout << "Smape error: " << smape_error << endl;

    return smape_error;
}


double compute_MSE(const vector<double>& res, const vector<double>& gt) {
    double sum_squared_errors = 0.0;
    int n = res.size();

    for (int i = 0; i < n; ++i) {
        double error = res[i] - gt[i];
        sum_squared_errors += error * error;
    }

    double mse = sum_squared_errors / n;
    cout << "MSE: " << mse << endl;
    return mse;
}

double compute_MSE_percentage(const vector<double>& res, const vector<double>& gt) {
    double sum_squared_errors = 0.0;
    double sum_squared_gt = 0.0;
    // int n = res.size();
    int n = 5344;

    for (int i = 0; i < n; ++i) {
        double error = res[i] - gt[i];
        sum_squared_errors += error * error;
        sum_squared_gt += gt[i] * gt[i];
    }

    if (sum_squared_gt == 0) {
        cerr << "Error: Sum of squared ground truth values is zero." << endl;
        return -1.0; // Return an invalid value to indicate error
    }

    double mse = sum_squared_errors / n;
    double mse_percentage = (mse / (sum_squared_gt / n)) * 100.0;

    cout << "MSE: " << mse << endl;
    cout << "MSE Percentage: " << mse_percentage << " %" << endl;

    return mse_percentage;
}


void run_ssgd_method(State &gs, int vertex, ofstream &smape_csv, ofstream &dist_csv) {
    // gs.res.clear();
    // cout << "Running SSGE method..." << endl;
    // gs.vtp_solver.query(vertex, gs.res);
    // cout << "SSGE method completed." << endl;
    // gs.res.resize(gs.nverts);

    gs.res.clear();
    cout << "Running Geotangle method..." << endl;
    gs.geotangle_solver.query(vertex, gs.res);
    cout << "Geotangle method completed." << endl;
    gs.res.resize(gs.nverts);

    // call the smape error function here
    double smape = smape_error(gs.res, gs.bob_ground_truth, gs);
    cout << "SMAPE: " << smape << endl;

    // Write the results to the csv file
    smape_csv << gs.mesh_name << "," << gs.nverts << "," << vertex << "," << smape << endl;

    for (int i = 0; i < gs.res.size(); i++) {
        dist_csv << i << "," << gs.res[i] << endl;
    }
}


// BACKUP
int main(int argc, char **argv) {
    State gs;

    // Valid vertices for the meshes
    // vector<int> vv_bob = {482, 1710, 2005, 3782, 4757};
    // vector<int> vv_bob = {482, 1710, 2005, 3782, 4757};
    vector<int> vv_bob = {482};


    vector<double> res_BOB_0 = vector<double>();
    vector<double> res_BOB_3 = vector<double>();
    vector<double> res_BOB_5 = vector<double>();
    
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

        string csv_BOB_0 = "../pymeshlab/Esperimento_1/data/bob_check/dist_BOB_0.csv";
        string csv_BOB_3 = "../pymeshlab/Esperimento_1/data/bob_check/dist_BOB_3.csv";
        string csv_BOB_5 = "../pymeshlab/Esperimento_1/data/bob_check/dist_BOB_5.csv";
        if (!read_ground_truth(csv_BOB_3, vertex, res_BOB_3)) {
            return 1;
        } else {
            cout << "Ground truth read successfully." << endl << endl;
        }

        // Verify that gs.bob_ground_truth and res_BOB_0 are not empty and print their sizes
        if (gs.bob_ground_truth.empty() || res_BOB_3.empty()) {
            cout << "Error: ground truth or res_BOB_3 is empty." << endl;
            return 1;
        } else {
            cout << "Ground truth size: " << gs.bob_ground_truth.size() << endl;
            cout << "Res_BOB_3 size: " << res_BOB_3.size() << endl;
        }

        // Compute the SMAPE error
        // double smape = smape_error(res_BOB_5, gs.bob_ground_truth, gs);
        double mse = compute_MSE(res_BOB_3, gs.bob_ground_truth);


    }

    return 0;
}

// BOB_3: no percentage -> 0.00814959
// BOB_5: no percentage -> 0.00828597


