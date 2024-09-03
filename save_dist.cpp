// *********************************************************
// BATCH FILE TO SAVE ALL THE DISTANCES IN .txt FILES
// *********************************************************

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
#include <filesystem>

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


const string OUTPUT_PATH = "../pymeshlab/Esperimento_2/prova";


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
    TrettnerSolver trettner_solver;
    FastMarchingSolver fast_mar_solver;
    HeatSolver heat_solver;
    GeotangleSolver geotangle_solver;
    EdgeSolver edge_solver;
    ExtendedSolver extended_solver;
    LanthierSolver lanthier_solver;

    // Timer
    double vtp_pre, tre_pre, fast_pre, heat_pre, geo_pre, edge_pre, ext_pre, lan_pre;

    // double preproc_time;
    double query_time;

    vector<int> sources;

    State() {
        // sources = {100}; // Default source vertex
        // Dragon source:
        // sources = {16};

        // Timer
        vtp_pre = 0.0, tre_pre = 0.0, fast_pre = 0.0, heat_pre = 0.0;
        geo_pre = 0.0, edge_pre = 0.0, ext_pre = 0.0, lan_pre = 0.0;

        // preproc_time = 0.0;
        query_time = 0.0;
    }
};

void Load_mesh(string filename, State &gs) {
    gs.m = DrawableTrimesh<>(filename.c_str());
    gs.nverts = gs.m.num_verts();
    gs.VV.resize(gs.nverts);

    gs.coords = extract_coords(gs.m);
    gs.tris = extract_tris(gs.m);

    for (auto i = 0; i < gs.nverts; i++)
        gs.VV[i] = gs.m.vert_ordered_verts_link(i);

    gs.m.normalize_bbox();
    gs.m.center_bbox();


    // FOR TRETTER
    // Save thi normalized mesh to a file .obj
    // Extract the parent path of the file
    string parent_path = fs::path(filename).parent_path().string();
    // Go back to the parent of the parent folder
    string grandparent_path = fs::path(parent_path).parent_path().string();
    // Construct the new path for the "normalized_for_trettner" folder
    string new_folder = grandparent_path + "/normalized_for_trettner";
    // Ensure the new folder exists
    fs::create_directories(new_folder);
    // Construct the new filename
    string new_filename = fs::path(filename).filename().string();
    string out_normalized_bb_mesh = new_folder + "/" + new_filename.substr(0, new_filename.size() - 4) + "_NORMALIZED.obj";
    // Save the normalized mesh
    gs.m.save(out_normalized_bb_mesh.c_str());
    // Give to Trettner the normalized mesh
    gs.trettner_solver = TrettnerSolver(out_normalized_bb_mesh);
    // gs.trettner_solver = TrettnerSolver(filename);
}

void init(GeodesicMethod &method, State &gs, const string &name) {
    cout << "----- Init method: " << name << " -----" << endl;
    cout << "----- Loading method..." << endl;
    method.load(&gs.m);

    cout << "----- Preprocessing method..." << endl;
    auto tic = std::chrono::steady_clock::now();
    method.preprocess();
    auto toc = std::chrono::steady_clock::now();
    double preprocess_time = std::chrono::duration_cast<std::chrono::duration<double>>(toc - tic).count();
    cout << endl << "$$$ PREPROCESS: " << preprocess_time << endl << endl;

    if (name == "VTP") gs.vtp_pre = preprocess_time;
    if (name == "Trettner") gs.tre_pre = preprocess_time;
    if (name == "Fast Marching") gs.fast_pre = preprocess_time;
    if (name == "Heat") gs.heat_pre = preprocess_time;
    if (name == "Geotangle") gs.geo_pre = preprocess_time;
    if (name == "Edge") gs.edge_pre = preprocess_time;
    if (name == "Extended") gs.ext_pre = preprocess_time;
    if (name == "Lanthier") gs.lan_pre = preprocess_time;

    
    // gs.preproc_time = preprocess_time;
}

void init_methods(State &gs) {
    init(gs.vtp_solver, gs, "VTP");
    init(gs.trettner_solver, gs, "Trettner");
    init(gs.fast_mar_solver, gs, "Fast Marching");
    init(gs.heat_solver, gs, "Heat");
    init(gs.geotangle_solver, gs, "Geotangle");
    init(gs.edge_solver, gs, "Edge");
    init(gs.lanthier_solver, gs, "Lanthier");
    // init(gs.extended_solver, gs, "Extended");
    
}

void write_results_to_file(const string& mesh_name, const string& method_name, int source_vertex,
                           const State& gs, const vector<double>& res) {

    fs::path file_path = fs::path(OUTPUT_PATH) / (mesh_name + "_" + method_name + "_" + to_string(source_vertex) + ".txt");
    ofstream outfile(file_path);

    double preproc_time = 0.0;
    // if (method_name == "VTP") preproc_time = gs.vtp_pre;
    // else if (method_name == "Trettner") preproc_time = gs.tre_pre;
    // else if (method_name == "Fast Marching") preproc_time = gs.fast_pre;
    if (method_name == "VTP") preproc_time = 0.0;
    else if (method_name == "Trettner") preproc_time = 0.0;
    else if (method_name == "Fast Marching") preproc_time = 0.0;
    else if (method_name == "Heat") preproc_time = gs.heat_pre;
    else if (method_name == "Geotangle") preproc_time = gs.geo_pre;
    else if (method_name == "Edge") preproc_time = gs.edge_pre;
    else if (method_name == "Extended") preproc_time = gs.ext_pre;
    else if (method_name == "Lanthier") preproc_time = gs.lan_pre;
    
    if (outfile.is_open()) {
        outfile << "# " << mesh_name << endl;
        outfile << "# Number of vertices: " << gs.nverts << ", Number of faces: " << gs.m.num_polys() << endl;
        outfile << "# Preprocessing time: " << preproc_time << endl;
        outfile << "# Query time: " << gs.query_time << endl;
        
        for (const auto& value : res) {
            outfile << value << endl;
        }
        
        outfile.close();
        cout << "Results written to " << file_path << endl;
    } else {
        cerr << "Unable to open file: " << file_path << endl;
    }
}

void run_ssgd_method(State &gs, int vertex, const string& mesh_name) {
    cout << "----- Running SSGD method -----" << endl;

auto run_method = [&](auto& solver, const string& method_name) {
        cout << "----- " << method_name << " -----" << endl;
        auto start = std::chrono::steady_clock::now();
        solver.query(vertex, gs.res);
        auto end = std::chrono::steady_clock::now();
        gs.query_time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
        
        write_results_to_file(mesh_name, method_name, vertex, gs, gs.res);
        gs.res.clear();
    };

    run_method(gs.vtp_solver, "VTP");
    run_method(gs.trettner_solver, "Trettner");
    run_method(gs.fast_mar_solver, "Fast Marching");
    run_method(gs.heat_solver, "Heat");
    run_method(gs.geotangle_solver, "Geotangle");
    run_method(gs.edge_solver, "Edge");
    run_method(gs.lanthier_solver, "Lanthier");
    // run_method(gs.extended_solver, "Extended");
}

int main(int argc, char **argv) {
    State gs;

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <folder_path>" << endl;
        return 1;
    }
    string folder_path = argv[1];

    // Create output directory if it doesn't exist
    fs::create_directories(OUTPUT_PATH);

    // vector<int> vv_sources = {100};
    // vector<int> vv_sources = {729, 1989, 2519, 93, 475};
    vector<int> vv_sources = {729};


    for(int vertex : vv_sources) {
        for (const auto &entry : fs::directory_iterator(folder_path)) {
            if (entry.path().extension() == ".obj") {
                string mesh_path = entry.path().string();
                string mesh_name = entry.path().stem().string();
                cout << "Processing mesh: " << mesh_path << endl;

                Load_mesh(mesh_path, gs);
            
                init_methods(gs);

                run_ssgd_method(gs, vertex, mesh_name);
            }
        } 
    }

    return 0;
}

/* Per usarlo: ./save_dist_app ../pymeshlab/Esperimento_Distribuzioni/data/bunny_500f
 * Dove nella cartella ci sono i file .obj su cui calcolare le distanze con i metodi SSGD   
 */