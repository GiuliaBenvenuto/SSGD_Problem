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
    DrawableTrimesh<> m;    // the input mesh
    uint nverts;            // number of vertices
    vector<uint> tris;      // triangle indices
    vector<double> coords;  // vertex coordinates
    vector<double> res;     // result
    ScalarField field;

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

    // vector<double> blub_ground_truth;
    // vector<double> bob_ground_truth;
    vector<double> spot_ground_truth;

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

        // blub_ground_truth = vector<double>();
        // bob_ground_truth = vector<double>();
        spot_ground_truth = vector<double>();

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

class MeshCache {
private:
    struct MeshData {
        DrawableTrimesh<> mesh;
        std::vector<uint> tris;
        std::vector<double> coords;
    };

    std::unordered_map<std::string, MeshData> cache;


public:
    bool isMeshLoaded(const std::string& path) {
        return cache.find(path) != cache.end();
    }

    MeshData& loadMesh(const std::string& path) {
        if (!isMeshLoaded(path)) {
            DrawableTrimesh<> m(path.c_str());
            cache[path] = {m, extract_tris(m), extract_coords(m)};
        }
        return cache[path];
    }
};


// Load meshes
void load_mesh(const string &filename, State &gs, MeshCache& cache) {
    if (!cache.isMeshLoaded(filename)) {
        cout << "LOADING MESH BECAUSE IT IS NOT IN CACHE" << endl;
        auto& meshData = cache.loadMesh(filename);
        gs.m = meshData.mesh;
        gs.nverts = gs.m.num_verts();
        gs.coords = meshData.coords;
        gs.tris = meshData.tris;
        gs.trettner_solver.mesh_path = filename;
        // Perform solver-specific initializations here
    } else {
        cout << "USING MESH FROM CACHE" << endl;
        auto& meshData = cache.loadMesh(filename);
        gs.m = meshData.mesh;
        gs.nverts = gs.m.num_verts();
        gs.coords = meshData.coords;
        gs.tris = meshData.tris;
        gs.trettner_solver.mesh_path = filename;
    } 
}


void init(GeodesicMethod &m, State &gs, const string &name) {
    cout << endl << "---------- Initializing method: " << name << " ----------" << endl;

    // Load
    gs.tic = chrono::steady_clock::now();
    m.load(&gs.m);
    gs.toc = chrono::steady_clock::now();
    double load_time = chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic).count();
    cout << "LOAD TIME: " << load_time << " milliseconds" << endl;

    if (name == "Heat") {   
        double time_scalar = 1;  // Adjust with a factor (e.g., 10.0)
        
        cout << "Time scalar: " << time_scalar << endl;
        gs.heat_solver.set_t(time_scalar);

    }

    // Preprocess
    gs.tic = chrono::steady_clock::now();
    m.preprocess();
    gs.toc = chrono::steady_clock::now();
    double preprocess_time = chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic).count();
    cout << "PREPROCESS TIME: " << preprocess_time << " milliseconds" << endl;

    // Update timings in the state
    if (name == "VTP") {
        gs.vtp_load = load_time;
        gs.vtp_preprocess = preprocess_time;
    } 
    else if (name == "Trettner") {
        gs.trettner_load = load_time;
        gs.trettner_preprocess = preprocess_time;
    }
    else if (name == "Fast Marching") {
        gs.fast_mar_load = load_time;
        gs.fast_mar_preprocess = preprocess_time;
    }
    else if (name == "Heat") {
        gs.heat_load = load_time;
        gs.heat_preprocess = preprocess_time;
    }
    else if (name == "Geotangle") {
        gs.geotangle_load = load_time;
        gs.geotangle_preprocess = preprocess_time;
    }
    else if (name == "Edge") {
        gs.edge_load = load_time;
        gs.edge_preprocess = preprocess_time;
    }
    else if (name == "Extended") {
        gs.extended_load = load_time;
        gs.extended_preprocess = preprocess_time;
    }
    else if (name == "Lanthier") {
        gs.lanthier_load = load_time;
        gs.lanthier_preprocess = preprocess_time;
    }
}

void init_methods(State &gs) {
    init(gs.vtp_solver,         gs,     "VTP");
    // init(gs.trettner_solver,    gs,     "Trettner");
    // init(gs.fast_mar_solver,    gs,     "Fast Marching");
    // init(gs.heat_solver,        gs,     "Heat");
    // init(gs.geotangle_solver,   gs,     "Geotangle");
    // init(gs.edge_solver,        gs,     "Edge");
    // init(gs.lanthier_solver,    gs,     "Lanthier");
    // init(gs.extended_solver,    gs,     "Extended");
}

// Function to calculate SMAPE between two vectors
double calculate_smape(const vector<double>& gt, const vector<double>& est) {
    cout << "GT size: " << gt.size() << ", EST size: " << est.size() << endl;
    if (gt.empty() || est.empty()) {
        cerr << "Ground truth or estimated distances are empty." << endl;
        return 0.0;
    }

    double smape = 0.0;
    int count = 0;

    // for (size_t i = 0; i < gt.size() && i < est.size(); ++i) {
    for (size_t i = 0; i < est.size(); ++i) {
        double denom = std::abs(gt[i]) + std::abs(est[i]);
        if (denom != 0) {
            smape += (std::abs(gt[i] - est[i]) / denom);
            ++count;
        }
    }

    if (count > 0) {
        smape = (smape / count) * 100.0;  // Convert to percentage
    }
    return smape;
}


void run_ssgd_method(State &state, int sourceVertexIndex, string type, vector<double> &gt, vector<double>& smape_errors) {
    // Initialize variables
    vector<double> distances = vector<double>(state.nverts, DBL_MAX);
    // vector<double> distances;

    ScalarField field;

    //vector<double> &ground_truth = gt;
    vector<double> ground_truth = gt;

    // check if ground truth is empty
    if (ground_truth.empty()) {
        cerr << "Ground truth is empty." << endl;
        return;
    } else {
        cout << "Ground truth is NOT empty." << endl;
    }

    auto log_time_and_calculate_smape = [&](auto &solver, const string &method) {
        auto start = chrono::high_resolution_clock::now();

        solver.query(sourceVertexIndex, distances);  // Only the query operation

        auto end = chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        cout << "QUERY TIME: " << elapsed.count() << " s" << endl;

        // Calculate SMAPE using ground_truth instead of gt
        if (method == "Lanthier") {
            // take only the first nverts values
            distances.resize(state.nverts);
            cout << "Distances size for LANTHIER after resize: " << distances.size() << endl;
        }
        double smape = calculate_smape(gt, distances);
        cout << method << " SMAPE: " << smape << "%" << endl;

        // Save SMAPE error for this method
        smape_errors.push_back(smape);
        distances.clear();
    };

    // VTP Solver
    cout << endl << "----- VTP -----" << endl;
    log_time_and_calculate_smape(state.vtp_solver, "VTP");

    // // Trettner Solver
    // cout << endl << "----- Trettner -----" << endl;
    // log_time_and_calculate_smape(state.trettner_solver, "Trettner");

    // // Fast Marching Solver
    // cout << endl << "----- Fast Marching Query -----" << endl;
    // log_time_and_calculate_smape(state.fast_mar_solver, "Fast Marching");

    // // Heat Solver
    // cout << endl << "----- Heat -----" << endl;
    // log_time_and_calculate_smape(state.heat_solver, "Heat");

    // // Geotangle Solver
    // cout << endl << "----- Geotangle -----" << endl;
    // log_time_and_calculate_smape(state.geotangle_solver, "Geotangle");

    // // Edge Solver
    // cout << endl << "----- Edge -----" << endl;
    // log_time_and_calculate_smape(state.edge_solver, "Edge");

    // // Lanthier Solver
    // cout << endl << "----- Lanthier -----" << endl;
    // log_time_and_calculate_smape(state.lanthier_solver, "Lanthier");

    // //Extended Solver
    // cout << endl << "----- Extended -----" << endl;
    // log_time_and_calculate_smape(state.extended_solver, "Extended");
}


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


// Assume other necessary headers and namespace declarations are here
int main(int argc, char **argv) {
    State gs;
    MeshCache cache;

    // Valid vertices for the meshes
    // vector<int> vv_blub = {663, 3958, 4662, 4715, 6694};
    // vector<int> vv_bob = {1710, 3782, 4757, 482, 2005};
    vector<int> vv_spot = {395, 2794, 283, 174, 1876}; 

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <folder_path>" << endl;
        return 1;
    }
    string folderPath = argv[1];

    for (int vertex : vv_spot) {
        cout << "------- Processing vertex: " << vertex << " --------" << endl;
        gs.spot_ground_truth.clear();

        // Prepare CSV file
        ofstream csvFile("../pymeshlab/Esperimento_1/data/smape/smape_spot_VTP_" + to_string(vertex) + ".csv");
        csvFile << "MeshName,NumVertices,SMAPE_vtp\n";

        // read a csv file
        string csv_gt = "../pymeshlab/Esperimento_1/data/gt/spot_gt_distances.csv";
        cout << "Reading ground truth from: " << csv_gt << endl;

        // Use the read_ground_truth function
        if (!read_ground_truth(csv_gt, vertex, gs.spot_ground_truth)) {
            cerr << "Failed to read ground truth for vertex " << vertex << endl;
            return 2;
        }

        for (const auto &entry : fs::directory_iterator(folderPath)) {
            if (entry.path().extension() == ".obj") {
                string meshPath = entry.path().string();
                gs.mesh_path = meshPath;
                gs.mesh_name = entry.path().filename().string();
                cout << "----- Mesh name: " << gs.mesh_name << " -----" << endl;
                string type = gs.mesh_name.substr(0, 4);

                vector<double> smape_errors;

                try {
                    load_mesh(meshPath, gs, cache);
                    init_methods(gs);
                    run_ssgd_method(gs, vertex, type, gs.spot_ground_truth, smape_errors);

                    // Write to CSV after processing each mesh
                    csvFile << gs.mesh_name << "," << gs.nverts;
                    for (const auto& smape : smape_errors) {
                        csvFile << "," << smape;
                    }
                    csvFile << "\n";
                    csvFile.flush(); // Ensure the line is written immediately

                    cout << "###### Results for " << gs.mesh_name << " written to CSV." << endl << endl;

                } catch (const std::exception& e) {
                    cerr << "Error processing mesh " << gs.mesh_name << ": " << e.what() << endl;
                    csvFile << gs.mesh_name << "," << gs.nverts << ",ERROR\n";
                    csvFile.flush();
                }
            }
        }
        csvFile.close(); // Ensure the CSV file is closed after all meshes for the vertex are processed
    }

    return 0;
}


// CON DUE FOLDER PATH
// int main(int argc, char **argv) {
//     State gs;
//     MeshCache cache;

//     // Valid vertices for the meshes
//     // vector<int> vv_blub = {663, 3958, 4662, 4715, 6694};
//     vector<int> vv_bob = {1710, 3782, 4757, 482, 2005};
//     // vector<int> vv_spot = {395, 2794, 283, 174, 1876}; 

//     if (argc < 3) {
//         cerr << "Usage: " << argv[0] << " <folder_path_1>" << " <folder_path_2>" << endl;
//         return 1;
//     }
//     string folderPath_1 = argv[1];
//     string folderPath_2 = argv[2];

//     // PRIMO GIRO SU BOB
//     for (int vertex : vv_bob) {
//         cout << "------- Processing vertex: " << vertex << " --------" << endl;
//         gs.bob_ground_truth.clear();

//         // Prepare CSV file
//         ofstream csvFile("../pymeshlab/Esperimento_1/data/smape/smape_bob_ext_" + to_string(vertex) + ".csv");
//         csvFile << "MeshName,NumVertices,SMAPE_Extended\n";

//         // read a csv file
//         string csv_gt = "../pymeshlab/Esperimento_1/data/gt/bob_gt_distances.csv";
//         cout << "Reading ground truth from: " << csv_gt << endl;

//         // read only the column of csv_gt that is names "vertex_" + to_string(vertex)
//         ifstream gt_csvFile(csv_gt);
//         string line;
//         getline(gt_csvFile, line); 
//         stringstream header(line);
//         string headerValue;
//         int columnIndex = -1;
//         int index = 0;

//         // Find the correct column for the vertex
//         while (getline(header, headerValue, ',')) {
//             if (headerValue == "vertex_" + to_string(vertex)) {
//                 columnIndex = index;
//                 break;
//             }
//             index++;
//         }

//         if (columnIndex == -1) {
//             cerr << "Vertex " << vertex << " not found in header." << endl;
//             return 2;
//         } else {
//             cout << "Vertex " << vertex << " found at column " << columnIndex << endl;
//         }

//         // Read the specified column from each line
//         while (getline(gt_csvFile, line)) {
//             stringstream lineStream(line);
//             string cell;
//             int currentColumn = 0;

//             while (getline(lineStream, cell, ',')) {
//                 if (currentColumn == columnIndex) {
//                     gs.bob_ground_truth.push_back(stod(cell));
//                     break;
//                 }
//                 currentColumn++;
//             }
//         }

//         if (gs.bob_ground_truth.empty()) {
//             cerr << "Ground truth is empty." << endl;
//             return 1;
//         } else {
//             cout << gs.bob_ground_truth.size() << " ground truth values read." << endl;
//         }

//         gt_csvFile.close();

//         for (const auto &entry : fs::directory_iterator(folderPath_1)) {
//             if (entry.path().extension() == ".obj") {
//                 string meshPath = entry.path().string();
//                 gs.mesh_path = meshPath;
//                 gs.mesh_name = entry.path().filename().string();
//                 cout << "----- Mesh name: " << gs.mesh_name << " -----" << endl;
//                 string type = gs.mesh_name.substr(0, 4);

//                 vector<double> smape_errors;

//                 try {
//                     load_mesh(meshPath, gs, cache);
//                     init_methods(gs);
//                     run_ssgd_method(gs, vertex, type, gs.bob_ground_truth, smape_errors);

//                     // Write to CSV after processing each mesh
//                     csvFile << gs.mesh_name << "," << gs.nverts;
//                     for (const auto& smape : smape_errors) {
//                         csvFile << "," << smape;
//                     }
//                     csvFile << "\n";
//                     csvFile.flush(); // Ensure the line is written immediately

//                     cout << "###### Results for " << gs.mesh_name << " written to CSV." << endl << endl;

//                 } catch (const std::exception& e) {
//                     cerr << "Error processing mesh " << gs.mesh_name << ": " << e.what() << endl;
//                     csvFile << gs.mesh_name << "," << gs.nverts << ",ERROR\n";
//                     csvFile.flush();
//                 }
//             }
//         }
//         csvFile.close(); // Ensure the CSV file is closed after all meshes for the vertex are processed
//     }



//     // SECONDO GIRO SU BLUB
//     vector<int> vv_blub = {663, 3958, 4662, 4715, 6694};

//     for (int vertex : vv_blub) {
//         cout << "------- Processing vertex: " << vertex << " --------" << endl;
//         gs.blub_ground_truth.clear();

//         // Prepare CSV file
//         ofstream csvFile("../pymeshlab/Esperimento_1/data/smape/smape_blub_ext_" + to_string(vertex) + ".csv");
//         csvFile << "MeshName,NumVertices,SMAPE_Extended\n";

//         // read a csv file
//         string csv_gt = "../pymeshlab/Esperimento_1/data/gt/blub_gt_distances.csv";
//         cout << "Reading ground truth from: " << csv_gt << endl;

//         // read only the column of csv_gt that is names "vertex_" + to_string(vertex)
//         ifstream gt_csvFile(csv_gt);
//         string line;
//         getline(gt_csvFile, line); 
//         stringstream header(line);
//         string headerValue;
//         int columnIndex = -1;
//         int index = 0;

//         // Find the correct column for the vertex
//         while (getline(header, headerValue, ',')) {
//             if (headerValue == "vertex_" + to_string(vertex)) {
//                 columnIndex = index;
//                 break;
//             }
//             index++;
//         }

//         if (columnIndex == -1) {
//             cerr << "Vertex " << vertex << " not found in header." << endl;
//             return 2;
//         } else {
//             cout << "Vertex " << vertex << " found at column " << columnIndex << endl;
//         }

//         // Read the specified column from each line
//         while (getline(gt_csvFile, line)) {
//             stringstream lineStream(line);
//             string cell;
//             int currentColumn = 0;

//             while (getline(lineStream, cell, ',')) {
//                 if (currentColumn == columnIndex) {
//                     gs.blub_ground_truth.push_back(stod(cell));
//                     break;
//                 }
//                 currentColumn++;
//             }
//         }

//         if (gs.blub_ground_truth.empty()) {
//             cerr << "Ground truth is empty." << endl;
//             return 1;
//         } else {
//             cout << gs.blub_ground_truth.size() << " ground truth values read." << endl;
//         }

//         gt_csvFile.close();

//         for (const auto &entry : fs::directory_iterator(folderPath_2)) {
//             if (entry.path().extension() == ".obj") {
//                 string meshPath = entry.path().string();
//                 gs.mesh_path = meshPath;
//                 gs.mesh_name = entry.path().filename().string();
//                 cout << "----- Mesh name: " << gs.mesh_name << " -----" << endl;
//                 string type = gs.mesh_name.substr(0, 4);

//                 vector<double> smape_errors;

//                 try {
//                     load_mesh(meshPath, gs, cache);
//                     init_methods(gs);
//                     run_ssgd_method(gs, vertex, type, gs.blub_ground_truth, smape_errors);

//                     // Write to CSV after processing each mesh
//                     csvFile << gs.mesh_name << "," << gs.nverts;
//                     for (const auto& smape : smape_errors) {
//                         csvFile << "," << smape;
//                     }
//                     csvFile << "\n";
//                     csvFile.flush(); // Ensure the line is written immediately

//                     cout << "###### Results for " << gs.mesh_name << " written to CSV." << endl << endl;

//                 } catch (const std::exception& e) {
//                     cerr << "Error processing mesh " << gs.mesh_name << ": " << e.what() << endl;
//                     csvFile << gs.mesh_name << "," << gs.nverts << ",ERROR\n";
//                     csvFile.flush();
//                 }
//             }
//         }
//         csvFile.close(); // Ensure the CSV file is closed after all meshes for the vertex are processed
//     }

//     return 0;
// }
