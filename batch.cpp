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
namespace fs = std::filesystem;
/*
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

// Load meshes
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
*/


// la load() e preprocess() della mesh la faccio solo una volta per ogni mesh
// Poi una volta che ho la mesh e il preprocess eseguito, posso fare query()

// prima cosa: calcola ground truth 
// Quando runno un metodo k sulla mesh m, calcolo già l'errore rispetto alla ground truth con SMAPE
// Così nei csv ho direttamente l'errore rispetto alla ground truth


// ------- VTP sulle 3 mesh su un solo vertice -------
/*
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

    vector<double> blub_gt_dist;
    vector<double> bob_gt_dist;
    vector<double> spot_gt_dist;

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

        blub_gt_dist = vector<double>();
        bob_gt_dist = vector<double>();
        spot_gt_dist = vector<double>();

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

// Load meshes
void load_mesh(const string &filename, State &gs) {
    gs.m = DrawableTrimesh<>(filename.c_str());
    gs.nverts = gs.m.num_verts();
    gs.coords = extract_coords(gs.m);
    gs.tris = extract_tris(gs.m);
    gs.trettner_solver.mesh_path = filename;
}


void init(GeodesicMethod &m, State &gs, const string &name) {
    cout << "---------- Initializing method: " << name << " ----------" << endl;

    // Load
    gs.tic = chrono::steady_clock::now();
    m.load(gs.coords, gs.tris);
    gs.toc = chrono::steady_clock::now();
    double load_time = chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic).count();
    cout << "LOAD TIME: " << load_time << " milliseconds" << endl;

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
    init(gs.trettner_solver,    gs,     "Trettner");
    init(gs.fast_mar_solver,    gs,     "Fast Marching");
    init(gs.heat_solver,        gs,     "Heat");
    init(gs.geotangle_solver,   gs,     "Geotangle");
    init(gs.edge_solver,        gs,     "Edge");
    init(gs.extended_solver,    gs,     "Extended");
    init(gs.lanthier_solver,    gs,     "Lanthier");
}


void run_ssgd_method(State &state, int sourceVertexIndex, ofstream &csvFile) {
    vector<double> distances;
    ScalarField field;

    auto log_time = [&](auto &solver, const string &method) {
        auto tic = chrono::high_resolution_clock::now();

        solver.query(sourceVertexIndex, distances, field);  // Only the query operation

        auto toc = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(toc - tic).count();
        double time_ms = duration / 1000.0;
        cout << "Query time: " << time_ms << " ms" << endl;

        csvFile << state.mesh_name << "," << state.nverts << "," << sourceVertexIndex << ","
                << method << ",Query," << time_ms << "\n";
    };

    auto log_distances = [&](const string &method) {
        for (size_t i = 0; i < distances.size(); ++i) {
            csvFile << state.mesh_name << "," << state.nverts << "," << sourceVertexIndex << ","
                    << method << ",,," << i << "," << distances[i] << "\n";
        }
    };

    // VTP Solver
    cout << endl << "----- VTP -----" << endl;
    log_time(state.vtp_solver, "VTP");
    log_distances("VTP");

    // Trettner Solver
    cout << endl << "----- Trettner -----" << endl;
    log_time(state.trettner_solver, "Trettner");
    log_distances("Trettner");

    // Fast Marching Solver
    cout << endl << "----- Fast Marching -----" << endl;
    log_time(state.fast_mar_solver, "Fast Marching");
    log_distances("Fast Marching");

    // Heat Solver
    cout << endl << "----- Heat -----" << endl;
    log_time(state.heat_solver, "Heat");
    log_distances("Heat");

    // Geotangle Solver
    cout << endl << "----- Geotangle -----" << endl;
    log_time(state.geotangle_solver, "Geotangle");
    log_distances("Geotangle");

    // Edge Solver
    cout << endl << "----- Edge -----" << endl;
    log_time(state.edge_solver, "Edge");
    log_distances("Edge");

    // Lanthier Solver
    cout << endl << "----- Lanthier -----" << endl;
    log_time(state.lanthier_solver, "Lanthier");
    log_distances("Lanthier");

    // Extended Solver
    cout << endl << "----- Extended -----" << endl;
    log_time(state.extended_solver, "Extended");
    log_distances("Extended");

}


void write_gt(const std::string& filename, const std::vector<double>& blub_distances, const std::vector<double>& bob_distances, const std::vector<double>& spot_distances) {
    std::ofstream file(filename);

    // Check if the file is open
    if (!file.is_open()) {
        std::cerr << "Failed to open file for writing." << std::endl;
        return;
    }

    // Write the header
    file << "vertices,blub_gt_vtp,bob_gt_vtp,spot_gt_vtp\n";

    // Assuming all vectors are of the same size
    // spot_distances.size() is the number of vertices of the biggest mesh
    for (size_t i = 0; i < spot_distances.size(); ++i) {
        file << i << "," << blub_distances[i] << "," << bob_distances[i] << "," << spot_distances[i] << "\n";
    }

    file.close();
}


int main(int argc, char **argv) {
    vector<int> vertices_idx = {1, 651, 1301, 1951, 2601, 3251, 3901, 4551, 5201, 5850};

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

    // VTP Ground Truth
    // from the folderPath take the biggest meshes
    string blub_gt = "blub_tri_subdiv_5_final.obj";
    string bob_gt = "bob_tri_subdiv_5_final.obj";
    string spot_gt = "spot_tri_subdiv_6_final.obj";

    // on this meshs execute only vtp and save result in blub_gt_dist, bob_gt_dist, spot_gt_dist:
    State gs;
    double time_ms;
    long duration;
    cout << "---------- Initializing Ground Truth ----------" << endl;

    // cout << "---------- Blub ----------" << endl;
    // gs.mesh_path = (folderPath + "/" + blub_gt);
    // gs.mesh_name = blub_gt;
    // load_mesh(gs.mesh_path, gs);
    // init(gs.vtp_solver, gs, "VTP");
    // auto start = chrono::high_resolution_clock::now();
    // gs.vtp_solver.query(sourceVertexIndex, gs.blub_gt_dist, gs.field);
    // auto stop = chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed = stop - start;
    // cout << "VTP computation: " << elapsed.count() << " s" << endl;


    // cout << "---------- Bob ----------" << endl;
    // gs.mesh_path = (folderPath + "/" + bob_gt);
    // gs.mesh_name = bob_gt;
    // load_mesh(gs.mesh_path, gs);
    // init(gs.vtp_solver, gs, "VTP");
    // auto start = chrono::high_resolution_clock::now();
    // gs.vtp_solver.query(sourceVertexIndex, gs.bob_gt_dist, gs.field);
    // auto stop = chrono::high_resolution_clock::now();
    // std::chrono::duration<double> elapsed = stop - start;
    // cout << "VTP computation: " << elapsed.count() << " s" << endl;


    cout << "---------- Spot ----------" << endl;
    gs.mesh_path = (folderPath + "/" + spot_gt);
    gs.mesh_name = spot_gt;
    load_mesh(gs.mesh_path, gs);
    init(gs.vtp_solver, gs, "VTP");
    auto start = chrono::high_resolution_clock::now();
    gs.vtp_solver.query(sourceVertexIndex, gs.spot_gt_dist, gs.field);
    auto stop = chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = stop - start;
    cout << "VTP computation: " << elapsed.count() << " s" << endl;

    // write_gt("pymeshlab/Esperimento_1/data/csv/ground_truth.csv", gs.blub_gt_dist, gs.bob_gt_dist, gs.spot_gt_dist);
    // write the csv file with ground truth
    write_gt("ground_truth.csv", gs.blub_gt_dist, gs.bob_gt_dist, gs.spot_gt_dist);



    // for (const auto &entry : fs::directory_iterator(folderPath)) {
    //     if (entry.path().extension() == ".obj") {
    //         string meshPath = entry.path().string();
    //         State gs;
    //         gs.mesh_path = meshPath;
    //         gs.mesh_name = entry.path().filename().string();
    //         load_mesh(meshPath, gs);
            
    //         // remove .obj from gs.mesh_name
    //         gs.mesh_name = gs.mesh_name.substr(0, gs.mesh_name.size() - 4);
    //         string resultsFilename = (csvDirPath / (gs.mesh_name + "_idx_1.csv")).string();
    //         ofstream csvFile(resultsFilename, ios::trunc); // Truncate mode to overwrite data for each run
    //         csvFile << "Mesh Name,Total Vertices,Source Vertex Index,Method,Step,Time (ms),Target Vertex Index,Distance\n";

    //         init_methods(gs);
    //         run_ssgd_method(gs, sourceVertexIndex, csvFile);
            
    //         csvFile.close();
    //     }
    // }

    return 0;
}


// la load() e preprocess() della mesh la faccio solo una volta per ogni mesh
// Poi una volta che ho la mesh e il preprocess eseguito, posso fare query()

// prima cosa: calcola ground truth 
// Quando runno un metodo k sulla mesh m, calcolo già l'errore rispetto alla ground truth con SMAPE
// Così nei csv ho direttamente l'errore rispetto alla ground truth
*/


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

    // Timer
    std::chrono::steady_clock::time_point tic;
    std::chrono::steady_clock::time_point toc;
    double vtp_load,        vtp_preprocess,         vtp_query;

    vector<double> blub_gt_dist;
    vector<double> bob_gt_dist;
    vector<double> spot_gt_dist;

    State() {
        // Timer
        vtp_load = vtp_preprocess = vtp_query = 0.0;

        res = vector<double>();

        blub_gt_dist = vector<double>();
        bob_gt_dist = vector<double>();
        spot_gt_dist = vector<double>();

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

    // Load
    gs.tic = chrono::steady_clock::now();
    m.load(&gs.m);
    gs.toc = chrono::steady_clock::now();
    double load_time = chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic).count();
    cout << "LOAD TIME: " << load_time << " milliseconds" << endl;

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
    } else {
        cerr << "Method not found!" << endl;
    }
}

void init_methods(State &gs) {
    init(gs.vtp_solver,         gs,     "VTP");
}


void run_ssgd_method(State &state, int sourceVertexIndex, ofstream &csvFile) {
    vector<double> distances;
    ScalarField field;

    auto log_time = [&](auto &solver, const string &method) {
        auto tic = chrono::high_resolution_clock::now();

        solver.query(sourceVertexIndex, distances);  // Only the query operation

        auto toc = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(toc - tic).count();
        double time_ms = duration / 1000.0;
        cout << "Query time: " << time_ms << " ms" << endl;

        csvFile << state.mesh_name << "," << state.nverts << "," << sourceVertexIndex << ","
                << method << ",Query," << time_ms << "\n";
    };

    auto log_distances = [&](const string &method) {
        for (size_t i = 0; i < distances.size(); ++i) {
            csvFile << state.mesh_name << "," << state.nverts << "," << sourceVertexIndex << ","
                    << method << ",,," << i << "," << distances[i] << "\n";
        }
    };

    // VTP Solver
    cout << endl << "----- VTP -----" << endl;
    log_time(state.vtp_solver, "VTP");
    log_distances("VTP");

}


// void write_gt(const std::string& filename, const std::vector<double>& blub_distances, const std::vector<double>& bob_distances, const std::vector<double>& spot_distances) {
//     std::ofstream file(filename);

//     // Check if the file is open
//     if (!file.is_open()) {
//         std::cerr << "Failed to open file for writing." << std::endl;
//         return;
//     }

//     // Write the header
//     file << "vertices,blub_gt_vtp,bob_gt_vtp,spot_gt_vtp\n";

//     // // spot_distances.size() is the number of vertices of the biggest mesh
//     // for (size_t i = 0; i < spot_distances.size(); ++i) {
//     //     file << i << "," << blub_distances[i] << "," << bob_distances[i] << "," << spot_distances[i] << "\n";
//     // }
//     for (size_t i = 0; i < blub_distances.size(); ++i) {
//         file << i << "," << blub_distances[i] << "," << bob_distances[i] << "," << spot_distances[i] << "\n";
//     }

//     file.close();
// }


// void write_gt(const std::string& filename, const std::vector<double>& blub_distances, const std::vector<double>& bob_distances, const std::vector<double>& spot_distances) {    
//     std::ofstream file(filename);

//     if (!file.is_open()) {
//         std::cerr << "Failed to open file for writing." << std::endl;
//         return;
//     }

//     file << "vertices,blub_gt_vtp,bob_gt_vtp,spot_gt_vtp\n";
//     size_t max_size = std::max({blub_distances.size(), bob_distances.size(), spot_distances.size()});

//     for (size_t i = 0; i < max_size; ++i) {
//         double blub_dist = i < blub_distances.size() ? blub_distances[i] : 0.0;
//         double bob_dist = i < bob_distances.size() ? bob_distances[i] : 0.0;
//         double spot_dist = i < spot_distances.size() ? spot_distances[i] : 0.0;
//         file << i << "," << blub_dist << "," << bob_dist << "," << spot_dist << "\n";
//     }
//     file.close();
// }

void write_gt(const std::string& filename, const std::vector<double>& blub_distances, const std::vector<double>& bob_distances, const std::vector<double>& spot_distances) {    
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Failed to open file for writing." << std::endl;
        return;
    }

    file << "vertices,blub_gt_vtp,bob_gt_vtp,spot_gt_vtp\n";
    size_t max_size = std::max({blub_distances.size(), bob_distances.size(), spot_distances.size()});
    cout << "max_size: " << max_size << endl;

    bool has_negative = false;
    for (size_t i = 0; i < max_size; ++i) {
        double blub_dist = i < blub_distances.size() ? blub_distances[i] : -1;
        double bob_dist = i < bob_distances.size() ? bob_distances[i] : -1;
        double spot_dist = i < spot_distances.size() ? spot_distances[i] : -1;

        if (blub_dist < 0 || bob_dist < 0 || spot_dist < 0) {
            std::cerr << "Negative value detected at index " << i << std::endl;
            has_negative = true;
        }

        file << i << "," << blub_dist << "," << bob_dist << "," << spot_dist << "\n";
    }
    file.close();

    if (has_negative) {
        std::cerr << "There are negative values in the data." << std::endl;
    }
}




# include <chrono>
using namespace std::chrono;
int main(int argc, char **argv) {
    // vector<int> vertices_idx = {1, 651, 1301, 1951, 2601};
    vector<int> vertices_idx = {1, 651};

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <folder_path>" << endl;
        return 1;
    }
    string folderPath = argv[1];

    // Create a new directory for CSV files
    fs::path currentPath = fs::current_path();
    fs::path csvDirPath = currentPath / "csv_results";
    if (!fs::exists(csvDirPath)) {
        fs::create_directory(csvDirPath);
        cout << "Created directory: " << csvDirPath << endl;
    } else {
        cout << "Directory already exists: " << csvDirPath << endl;
    }

    // Mesh file names
    // string blub_gt = "blub_tri_subdiv_5_final.obj";
    // string bob_gt = "bob_tri_subdiv_5_final.obj";
    // string spot_gt = "spot_tri_subdiv_6_final.obj";

    string blub_gt = "blub_tri_subdiv_2_final.obj";
    string bob_gt = "bob_tri_subdiv_2_final.obj";
    string spot_gt = "spot_tri_subdiv_2_final.obj";

    State gs;
    cout << "---------- Initializing Ground Truth ----------" << endl;

    // Iterate over each vertex index
    for (int sourceVertexIndex : vertices_idx) {
        // Process each mesh for the current source vertex
        vector<string> meshes = {blub_gt, bob_gt, spot_gt};
        vector<vector<double>> results(3); // To store distances for blub, bob, spot

        for (size_t i = 0; i < meshes.size(); ++i) {
            const auto& mesh_file = meshes[i];
            gs.mesh_path = folderPath + "/" + mesh_file;
            gs.mesh_name = mesh_file;
            load_mesh(gs.mesh_path, gs);
            init(gs.vtp_solver, gs, "VTP");

            results[i].resize(gs.nverts);
            cout << "---> mesh: " << gs.mesh_name << "size: " << results[i].size() << endl;

            steady_clock::time_point start = steady_clock::now();
            gs.vtp_solver.query(sourceVertexIndex, results[i]);
            steady_clock::time_point end = steady_clock::now();
            // time in seconds
            double query_time = duration_cast<seconds>(end - start).count();

            cout << "MESH NAME: " << gs.mesh_name << ", VERTEX INDEX: " << sourceVertexIndex
                 << ", QUERY TIME: " << query_time << " ms" << endl << endl;
        }

        // Filename based on the vertex index
        string csvFileName = csvDirPath.string() + "/results_vertex_" + to_string(sourceVertexIndex) + ".csv";

        cout << "results[0].size(): " << results[0].size() << endl;
        cout << "results[1].size(): " << results[1].size() << endl;
        cout << "results[2].size(): " << results[2].size() << endl;

        write_gt(csvFileName, results[0], results[1], results[2]);
        cout << "Written distances for vertex index " << sourceVertexIndex << " to " << csvFileName << endl;
    }

    return 0;
}



/*
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

    vector<double> blub_ground_truth, bob_ground_truth, spot_ground_truth;

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

        blub_ground_truth = vector<double>();
        bob_ground_truth = vector<double>();
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

// Load meshes
void load_mesh(const string &filename, State &gs) {
    gs.m = DrawableTrimesh<>(filename.c_str());
    gs.nverts = gs.m.num_verts();
    gs.coords = extract_coords(gs.m);
    gs.tris = extract_tris(gs.m);
    gs.trettner_solver.mesh_path = filename;
}


void init(GeodesicMethod &m, State &gs, const string &name) {
    cout << endl << "---------- Initializing method: " << name << " ----------" << endl;

    // Load
    gs.tic = chrono::steady_clock::now();
    m.load(&gs.m);
    gs.toc = chrono::steady_clock::now();
    double load_time = chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic).count();
    cout << "LOAD TIME: " << load_time << " milliseconds" << endl;

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
    init(gs.trettner_solver,    gs,     "Trettner");
    init(gs.fast_mar_solver,    gs,     "Fast Marching");
    init(gs.heat_solver,        gs,     "Heat");
    init(gs.geotangle_solver,   gs,     "Geotangle");
    init(gs.edge_solver,        gs,     "Edge");
    // init(gs.extended_solver,    gs,     "Extended");
    init(gs.lanthier_solver,    gs,     "Lanthier");
}

// Funzione per calcolare SMAPE tra due vettori
// Valore alto di SMAPE indica che le due serie sono molto diverse
double calculate_smape(const vector<double>& gt, const vector<double>& est) {
    // cout << "GT size: " << gt.size() << ", EST size: " << est.size() << endl;
    // if (gt.empty() || est.empty()) {
    //     cerr << "Ground truth or estimated distances are empty." << endl;
    //     return 0.0;
    // }

    double smape = 0.0;
    int count = 0;

    for (size_t i = 0; i < gt.size() && i < est.size(); ++i) {
        double denom = std::abs(gt[i]) + std::abs(est[i]);
        if (denom != 0) {
            smape += std::abs(gt[i] - est[i]) / denom;
            ++count;
        }
    }

    if (count > 0) {
        smape = (smape / count) * 100.0;  // Convert to percentage
    }
    return smape;
}


void run_ssgd_method(State &state, int sourceVertexIndex, string type, vector<double> &gt) {
    vector<double> distances;
    ScalarField field;
    vector<double> &ground_truth = gt;


    // check if ground truth is empty
    if (ground_truth.empty()) {
        cerr << "Ground truth is empty." << endl;
        return;
    } else {
        cout << "Ground truth is NOT empty." << endl;
    }

    auto log_time_and_calculate_smape = [&](auto &solver, const string &method) {
        auto tic = chrono::high_resolution_clock::now();

        solver.query(sourceVertexIndex, distances);  // Only the query operation

        auto toc = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(toc - tic).count();
        double time_ms = duration / 1000.0;
        cout << "Query time: " << time_ms << " ms" << endl;

        // Calcolo di SMAPE using ground_truth instead of gt
        double smape = calculate_smape(ground_truth, distances);
        cout << method << " SMAPE: " << smape << "%" << endl;
    };

    // VTP Solver
    cout << endl << "----- VTP -----" << endl;
    log_time_and_calculate_smape(state.vtp_solver, "VTP");

    // Trettner Solver
    cout << endl << "----- Trettner -----" << endl;
    log_time_and_calculate_smape(state.trettner_solver, "Trettner");

    // Fast Marching Solver
    cout << endl << "----- Fast Marching Query -----" << endl;
    log_time_and_calculate_smape(state.fast_mar_solver, "Fast Marching");

    // Heat Solver
    cout << endl << "----- Heat -----" << endl;
    log_time_and_calculate_smape(state.heat_solver, "Heat");

    // Geotangle Solver
    cout << endl << "----- Geotangle -----" << endl;
    log_time_and_calculate_smape(state.geotangle_solver, "Geotangle");

    // Edge Solver
    cout << endl << "----- Edge -----" << endl;
    log_time_and_calculate_smape(state.edge_solver, "Edge");

    // Lanthier Solver
    cout << endl << "----- Lanthier -----" << endl;
    log_time_and_calculate_smape(state.lanthier_solver, "Lanthier");

    // Extended Solver
    // cout << endl << "----- Extended -----" << endl;
    // log_time_and_calculate_smape(state.extended_solver, "Extended");

}



// Assume other necessary headers and namespace declarations are here
# include <chrono>
using namespace std::chrono;
int main(int argc, char **argv) {
    State gs;
    // vector<int> vertices_idx = {1, 651, 1301, 1951, 2601, 3251, 3901, 4551, 5201, 5850};
    vector<int> vertices_idx = {1, 651, 1301, 1951, 2601};

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <folder_path>" << endl;
        return 1;
    }
    string folderPath = argv[1];

    // Create a new directory for CSV files
    fs::path currentPath = fs::current_path();
    fs::path csvDirPath = currentPath / "errors";
    if (!fs::exists(csvDirPath)) {
        fs::create_directory(csvDirPath);
        cout << "Created directory: " << csvDirPath << endl;
    } else {
        cout << "Directory already exists: " << csvDirPath << endl;
    }

    // ----- ORA VERTEX = 1 -----
    int vertex = vertices_idx[0];

    // read a csv file
    string csv_vtp_gt = "../pymeshlab/Esperimento_1/data/csv_vtp_gt/results_vertex_" + to_string(vertex) + ".csv";
    ifstream csvFile(csv_vtp_gt);
    // vector<double> blub_ground_truth, bob_ground_truth, spot_ground_truth;
    string line;
    getline(csvFile, line); // skip header

    while (getline(csvFile, line)) {
        std::stringstream ss(line);
        std::string cell;
        std::vector<double> temp;

        while (std::getline(ss, cell, ',')) {
            temp.push_back(std::stod(cell));
        }
        
        if (temp.size() > 3) {
            gs.blub_ground_truth.push_back(temp[1]); 
            gs.bob_ground_truth.push_back(temp[2]);      
            gs.spot_ground_truth.push_back(temp[3]);      
        }
    }

    // check if ground truth is empty
    if (gs.blub_ground_truth.empty() || gs.bob_ground_truth.empty() || gs.spot_ground_truth.empty()) {
        cerr << "Ground truth is empty." << endl;
        return 1;
    }


    // for each file .obj in the folderPath
    for (const auto &entry : fs::directory_iterator(folderPath)) {
        if (entry.path().extension() == ".obj") {
            string meshPath = entry.path().string();
            gs.mesh_path = meshPath;
            gs.mesh_name = entry.path().filename().string();
            string type = gs.mesh_name.substr(0, 4);


            load_mesh(meshPath, gs);
            init_methods(gs);

            if (type == "blub") {
                cout << "Blub size of ground truth: " << gs.blub_ground_truth.size() << endl;
                run_ssgd_method(gs, vertex, type, gs.blub_ground_truth);
            } else if (type == "bob_") {
                cout << "Bob size of ground truth: " << gs.bob_ground_truth.size() << endl;
                run_ssgd_method(gs, vertex, type, gs.bob_ground_truth);
            } else if (type == "spot") {
                cout << "Spot size of ground truth: " << gs.spot_ground_truth.size() << endl;
                run_ssgd_method(gs, vertex, type, gs.spot_ground_truth);
            } else {
                cerr << "Invalid type: " << type << endl;
                return 1;
            }
        }
    }

    return 0;
}
*/