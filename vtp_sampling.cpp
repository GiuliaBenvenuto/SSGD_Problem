/* CPP Program to find 5 valid vertices for which VTP works correctly */

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
#include <unordered_map>
#include <string>
#include <memory>

// Compute SSGD
#include "solving_ssgd.h"

using namespace std;
using namespace cinolib;
using namespace std::chrono;
namespace fs = std::filesystem;


struct State {
    DrawableTrimesh<> m;    // the input mesh
    uint nverts;            // number of vertices
    vector<uint> tris;      // triangle indices
    vector<double> coords;  // vertex coordinates
    vector<double> res;     // result
    ScalarField field;
    vector<int> sources;

    string mesh_path;
    string mesh_name;

    // Solver
    VTPSolver           vtp_solver;

    // Timer
    std::chrono::steady_clock::time_point tic;
    std::chrono::steady_clock::time_point toc;
    double vtp_load,        vtp_preprocess,         vtp_query;

    State() {
        // Timer
        vtp_load = vtp_preprocess = vtp_query = 0.0;

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
    gs.m.normalize_bbox();
    gs.m.center_bbox();
}

void init(GeodesicMethod &m, State &gs, const string &name) {
    // cout << endl << "---------- Initializing method: " << name << " ----------" << endl;
    // Load
    gs.tic = chrono::steady_clock::now();
    m.load(&gs.m);
    gs.toc = chrono::steady_clock::now();
    double load_time = chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic).count();
    //cout << "LOAD TIME: " << load_time << " milliseconds" << endl;

    // Preprocess
    gs.tic = chrono::steady_clock::now();
    m.preprocess();
    gs.toc = chrono::steady_clock::now();
    double preprocess_time = chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic).count();
    //cout << "PREPROCESS TIME: " << preprocess_time << " milliseconds" << endl;
}


void init_methods(State &gs) {
    init(gs.vtp_solver, gs, "VTP");
}


bool run_ssgd_method(State &state, int sourceVertexIndex) {
    ScalarField field;

    auto log_time_and_calculate_smape = [&](auto &solver, const string &method) {
        auto start = chrono::high_resolution_clock::now();
        cout << "QUERYING...";
        solver.query(sourceVertexIndex, state.res);  // Only the query operation

        auto end = chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        cout << "QUERY TIME: " << elapsed.count() << " s" << endl;

        for (double distance : state.res) {
            if (distance > 1000) {
                cout << "!!!!!!! INVALID DISTANCE: " << distance << endl;
                return false;
            }
        } 
        return true;
    };

    // VTP Solver
    cout << endl << "----- VTP -----" << endl;
    return log_time_and_calculate_smape(state.vtp_solver, "VTP");
}


int main(int argc, char **argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <folder_path>" << endl;
        return 1;
    }
    string folderPath = argv[1];
    State gs;
    vector<int> valid_vertices;

    // Load the reference mesh to sample vertices from
    // string refMeshPath = folderPath + "/blub_tri_final.obj";
    // string refMeshPath = folderPath + "/bob_tri_final.obj";
    // string refMeshPath = folderPath + "/M0_reordered.obj";
    string refMeshPath = folderPath + "/thai_ord_0.obj";
    
    load_mesh(refMeshPath, gs);
    cout << "Number of vertices in reference mesh: " << gs.nverts << endl;

    // ---- TAKE 50 RANDOM INDICES ----
    vector<int> indices(gs.nverts);
    iota(indices.begin(), indices.end(), 0);
    random_device rd;
    mt19937 g(rd());
    shuffle(indices.begin(), indices.end(), g);
    vector<int> selected_indices(indices.begin(), indices.begin() + 50);

    // VERTICI CHE VOGLIO TESTARE
    // selected_indices = {16, 167, 163, 194, 66};

    cout << "Random indices: ";
    for (int idx : selected_indices) {
        cout << idx << " ";
    }
    cout << endl;

    

    // Collect all .obj files in the folder
    vector<string> objFiles;
    for (const auto &entry : fs::directory_iterator(folderPath)) {
        if (entry.path().extension() == ".obj") {
            objFiles.push_back(entry.path().string());
        }
    }
    cout << "======= OBJ Files =======" << endl;
    for (const string &objPath : objFiles) {
        cout << objPath << endl;
    }

    // For each vertex selected
    for (int idx : selected_indices) {
        cout << "********* Vertex: " << idx << "*********" << endl;
        bool valid = true;

        // For each .obj file
        for (const string &objPath : objFiles) {
            cout << "Processing OBJ: " << objPath << endl;
            
            // Load the mesh
            load_mesh(objPath, gs);

            // Initialize the methods
            init_methods(gs);

            if (!run_ssgd_method(gs, idx)) {
                valid = false;
                cout << "Vertex " << idx << " is not valid for mesh: " << objPath << endl;
                break;
            } 
        }
        
        if (valid) {
            valid_vertices.push_back(idx);
            cout << "Vertex " << idx << " is valid for all meshes." << endl << endl;
            if (valid_vertices.size() == 5) {
                cout << "********* FOUND 5 VALID VERTICES *********" << endl;
                cout << "Valid vertices: ";
                for (int v : valid_vertices) {
                    cout << v << " ";
                }
                cout << endl;
                break;
            }
        } else {
            cout << "Vertex " << idx << " is not valid for all meshes." << endl;
        }
    }

    return 0;
}