// *******************************************
// QUESTO E' QUELLO PER CALCOLARE SMAPE GIUSTO 
// *******************************************

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

  double heat_time;
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

  vector<int> sources;

  State() {
    // sources = {100}; // Default source vertex
    heat_time = 1;
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

    // QUIIIII
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
}


void init(GeodesicMethod &method, State &gs, const string &name) {
    cout << "----- Init method: " << name << " -----" << endl;
    method.load(&gs.m);
    if (name == "Heat") {       
    gs.heat_solver.set_t(gs.heat_time);
    }
    method.preprocess();
}

void init_methods(State &gs) {
    // init(gs.vtp_solver, gs, "VTP");
    // init(gs.trettner_solver, gs, "Trettner");
    // init(gs.fast_mar_solver, gs, "Fast Marching");
    // init(gs.heat_solver, gs, "Heat");
    init(gs.geotangle_solver, gs, "Geotangle");
    // init(gs.edge_solver, gs, "Edge");
    init(gs.extended_solver, gs, "Extended");
    // init(gs.lanthier_solver, gs, "Lanthier");
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

    double smape = 0.0;
    int count = 0;

    for (int i = 0; i < res.size(); i++) {
        double num = abs(res[i] - gt[i]);
        double den = (abs(res[i]) + abs(gt[i])) / 2;
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

    return (smape / count) * 100.0;
}


void run_ssgd_method(State &gs, int vertex, vector<double> &ground_truth, vector<double> &smape_errors) {
    double smape = 0.0;

    cout << "----- Running SSGD method -----" << endl;

    // cout << "----- VTP ----- " << endl;
    // gs.vtp_solver.query(vertex, gs.res);
    // // Output results
    // // for (int i = 0; i < gs.res.size(); i++) {
    // //     cout << i << "," << gs.res[i] << endl;
    // // }
    // smape = smape_error(gs.res, ground_truth, gs);
    // smape_errors.push_back(smape);
    // gs.res.clear();
    // smape = 0.0;


    // cout << "----- TRETTNER -----" << endl;
    // gs.trettner_solver.query(vertex, gs.res);
    // // for (int i = 0; i < gs.res.size(); i++) {
    // //     cout << i << "," << gs.res[i] << endl;
    // // }
    // smape = smape_error(gs.res, ground_truth, gs);
    // smape_errors.push_back(smape);
    // gs.res.clear();
    // smape = 0.0;


    // cout << "----- FAST MARCHING -----" << endl;
    // gs.fast_mar_solver.query(vertex, gs.res);
    // smape = smape_error(gs.res, ground_truth, gs);
    // smape_errors.push_back(smape);
    // gs.res.clear();
    // smape = 0.0;


    // cout << "----- HEAT -----" << endl;
    // gs.heat_solver.query(vertex, gs.res);
    // smape = smape_error(gs.res, ground_truth, gs);
    // smape_errors.push_back(smape);
    // gs.res.clear();
    // smape = 0.0;


    cout << "----- GEOTANGLE -----" << endl;
    gs.geotangle_solver.query(vertex, gs.res);
    smape = smape_error(gs.res, ground_truth, gs);
    smape_errors.push_back(smape);
    gs.res.clear();
    smape = 0.0;


    // cout << "----- EDGE -----" << endl;
    // gs.edge_solver.query(vertex, gs.res);
    // smape = smape_error(gs.res, ground_truth, gs);
    // smape_errors.push_back(smape);
    // gs.res.clear();
    // smape = 0.0;

    
    cout << "----- EXTENDED -----" << endl;
    gs.extended_solver.query(vertex, gs.res);
    smape = smape_error(gs.res, ground_truth, gs);
    smape_errors.push_back(smape);
    gs.res.clear();
    smape = 0.0;


    // cout << "----- LANTHIER -----" << endl;
    // gs.lanthier_solver.query(vertex, gs.res);
    // gs.res.resize(gs.nverts);
    // smape = smape_error(gs.res, ground_truth, gs);
    // smape_errors.push_back(smape);
    // gs.res.clear();
    // smape = 0.0;

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


/*
// ----- MAIN FOR THE QUERY ------
int main(int argc, char **argv) {
    State gs;

    // TO RUN: ./bunny_app ../pymeshlab/Esperimento_1/data/bunny_prova
    if (argc < 2) {
    cerr << "Usage: " << argv[0] << " <folder_path>" << endl;
    return 1;
    }
    string folder_path = argv[1];

    // Read the ground truth from the CSV file we just created
    // std::filesystem::path csv_path = "../pymeshlab/Esperimento_1/data/gt/bunny_gt_distances.csv";
    // std::filesystem::path csv_path = "../pymeshlab/Esperimento_1/data/gt/bob500f_gt_distances.csv";
    // std::filesystem::path csv_path = "../pymeshlab/Esperimento_1/data/gt/Trettner_gt/bunny_gt_trettner.csv";
    // std::filesystem::path csv_path = "../pymeshlab/Esperimento_1_Thai/ground_truth/thai_decimated_gt.csv";
    // std::filesystem::path csv_path = "../pymeshlab/Esperimento_1/data/gt/spot500f_gt_distances.csv";
    std::filesystem::path csv_path = "../pymeshlab/Esperimento_1_Drago/Dragon_with_decimation/ground_truth/dragon_decimated_gt.csv";

    // vector<int> vv_bob = {100};
    // vector<int> vv_spot = {100};
    // vector<int> vv_spot500 = {100};
    // vector<int> vv_spot = {174};
    vector<int> vv_dragon = {16};
    // vector<int> vv_dragon = {16, 66, 163, 167, 194};
    // vector<int> vv_spot = {174, 283, 395, 1876, 2794};
    // vector<int> vv_bob = {482, 1710, 2005, 3782, 4757};
    // vector<int> vv_thai = {2, 48, 49, 50, 56};
    

    vector<double> ground_truth;

    for(int vertex : vv_dragon) {

      ground_truth.clear();

      // Prepare CSV file
      // ofstream csvFile("../pymeshlab/Esperimento_1/data/simplified/smape_BUNNY_TRETTNER_" + to_string(vertex) + ".csv");
      // ofstream csvFile("../pymeshlab/Esperimento_Drago/Dragon_with_decimation/SMAPE/smape_DRAGON_extended_" + to_string(vertex) + ".csv");
      // ofstream csvFile("../pymeshlab/Esperimento_1/data/SMAPE_RESULTS_BOB/smape_spot_VTP_" + to_string(vertex) + ".csv");
      ofstream csvFile("../pymeshlab/Esperimento_1_Drago/Dragon_with_decimation/SMAPE/smape_drago_HEAT_" + to_string(vertex) + ".csv");


      // csvFile << "MeshName,NumVertices,SMAPE_VTP,SMAPE_Trettner,SMAPE_FastMarching,SMAPE_Heat,SMAPE_Geotangle,SMAPE_Edge,SMAPE_Lanthier" << endl;
      csvFile << "MeshName,NumVertices,SMAPE_Heat" << endl;

      if (!read_ground_truth(csv_path.string(), vertex, ground_truth)) {
        cerr << "Failed to read ground truth from " << csv_path << endl;
        return 1;
      } else {
          cout << "Ground truth read successfully." << endl;
      }

      for (const auto &entry : fs::directory_iterator(folder_path)) {
            if (entry.path().extension() == ".obj") {
              string mesh_path = entry.path().string();
              cout << "Processing mesh: " << mesh_path << endl;

              Load_mesh(mesh_path, gs);
              // init_vtp(gs);
              init_methods(gs);

              vector<double> smape_errors;
              run_ssgd_method(gs, vertex, ground_truth, smape_errors);

              // Write the csv_file
              csvFile << entry.path().filename().string() << "," << gs.nverts << ",";
              for (double smape : smape_errors) {
                csvFile << smape << ",";
              }
              csvFile << "\n";
              csvFile.flush();
            }
      }
      csvFile.close();



    } // end for each vertex


    return 0;
}

*/


int main(int argc, char **argv) {
    State gs;

    string folder_path = "";
    std::filesystem::path csv_path = "";

    // ------- DRAGO PER DECIMAZIONE -------
    folder_path = "../pymeshlab/Esperimento_1_Drago/Dragon_with_decimation/dragon_ordered";

    // Read the ground truth from the CSV file we just created
    // std::filesystem::path csv_path = "../pymeshlab/Esperimento_1/data/gt/bunny_gt_distances.csv";
    // std::filesystem::path csv_path = "../pymeshlab/Esperimento_1/data/gt/bob500f_gt_distances.csv";
    // std::filesystem::path csv_path = "../pymeshlab/Esperimento_1/data/gt/Trettner_gt/bunny_gt_trettner.csv";
    // std::filesystem::path csv_path = "../pymeshlab/Esperimento_1_Thai/ground_truth/thai_decimated_gt.csv";
    // std::filesystem::path csv_path = "../pymeshlab/Esperimento_1/data/gt/spot500f_gt_distances.csv";
    csv_path = "../pymeshlab/Esperimento_1_Drago/Dragon_with_decimation/ground_truth/dragon_decimated_gt.csv";

    // vector<int> vv_bob = {100};
    // vector<int> vv_spot = {100};
    // vector<int> vv_spot500 = {100};
    // vector<int> vv_spot = {174};
    // vector<int> vv_dragon = {16};
    vector<int> vv_dragon = {16, 66, 163, 167, 194};
    // vector<int> vv_spot = {174, 283, 395, 1876, 2794};
    // vector<int> vv_bob = {482, 1710, 2005, 3782, 4757};
    // vector<int> vv_thai = {2, 48, 49, 50, 56};
    

    vector<double> ground_truth;

    for(int vertex : vv_dragon) {

      ground_truth.clear();
      ofstream csvFile("../pymeshlab/Esperimento_1_Drago/Dragon_with_decimation/SMAPE_EXT_GEO/smape_drago_" + to_string(vertex) + ".csv");


      // csvFile << "MeshName,NumVertices,SMAPE_VTP,SMAPE_Trettner,SMAPE_FastMarching,SMAPE_Heat,SMAPE_Geotangle,SMAPE_Edge,SMAPE_Lanthier" << endl;
      csvFile << "MeshName,NumVertices,SMAPE_Geotangle,SMAPE_Extended" << endl;

      if (!read_ground_truth(csv_path.string(), vertex, ground_truth)) {
        cerr << "Failed to read ground truth from " << csv_path << endl;
        return 1;
      } else {
          cout << "Ground truth read successfully." << endl;
      }

      for (const auto &entry : fs::directory_iterator(folder_path)) {
            if (entry.path().extension() == ".obj") {
              string mesh_path = entry.path().string();
              cout << "Processing mesh: " << mesh_path << endl;

              Load_mesh(mesh_path, gs);
              // init_vtp(gs);
              init_methods(gs);

              vector<double> smape_errors;
              run_ssgd_method(gs, vertex, ground_truth, smape_errors);

              // Write the csv_file
              csvFile << entry.path().filename().string() << "," << gs.nverts << ",";
              for (double smape : smape_errors) {
                csvFile << smape << ",";
              }
              csvFile << "\n";
              csvFile.flush();
            }
      }
      csvFile.close();



    } // end for each vertex


    // ------- THAI -------
    folder_path = "../pymeshlab/Esperimento_1_Thai/thai_ordered";
    csv_path = "../pymeshlab/Esperimento_1_Thai/ground_truth/thai_decimated_gt.csv";
    vector<int> vv_thai = {2, 48, 49, 50, 56};
    

    ground_truth.clear();

    for(int vertex : vv_thai) {

      ground_truth.clear();
      ofstream csvFile("../pymeshlab/Esperimento_1_Thai/SMAPE_EXT_GEO/smape_thai_" + to_string(vertex) + ".csv");


      // csvFile << "MeshName,NumVertices,SMAPE_VTP,SMAPE_Trettner,SMAPE_FastMarching,SMAPE_Heat,SMAPE_Geotangle,SMAPE_Edge,SMAPE_Lanthier" << endl;
      csvFile << "MeshName,NumVertices,SMAPE_Geotangle,SMAPE_Extended" << endl;

      if (!read_ground_truth(csv_path.string(), vertex, ground_truth)) {
        cerr << "Failed to read ground truth from " << csv_path << endl;
        return 1;
      } else {
          cout << "Ground truth read successfully." << endl;
      }

      for (const auto &entry : fs::directory_iterator(folder_path)) {
            if (entry.path().extension() == ".obj") {
              string mesh_path = entry.path().string();
              cout << "Processing mesh: " << mesh_path << endl;

              Load_mesh(mesh_path, gs);
              // init_vtp(gs);
              init_methods(gs);

              vector<double> smape_errors;
              run_ssgd_method(gs, vertex, ground_truth, smape_errors);

              // Write the csv_file
              csvFile << entry.path().filename().string() << "," << gs.nverts << ",";
              for (double smape : smape_errors) {
                csvFile << smape << ",";
              }
              csvFile << "\n";
              csvFile.flush();
            }
      }
      csvFile.close();



    } // end for each vertex


    return 0;
}


