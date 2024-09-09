
#include <Eigen/SparseCholesky>
#include <chrono>
#include <cinolib/color.h>
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
#include <future>
#include <imgui.h>
#include <thread>

// Compute SSGD
#include "solving_ssgd.h"

// Matlab
#include "MatlabDataArray.hpp"
#include "MatlabEngine.hpp"

using namespace std;
using namespace cinolib;
using namespace gcHeatWrapper;
using namespace matlab::engine;
namespace fs = std::filesystem;

//------ Global variables ------
// Fonts
ImFont *lato_bold = nullptr;
ImFont *lato_regular = nullptr;
ImFont *lato_bold_title = nullptr;
atomic<float> progress(0.0f);

//:::::::::::::::::::::::::::: GLOBAL VARIABLES (FOR
//: GUI):::::::::::::::::::::::::::::::
struct State {
  //-------- Program state --------
  bool MESH_IS_LOADED;
  // Input
  DrawableTrimesh<> m;     // the input mesh
  uint nverts;             // its #vertices
  vector<vector<uint>> VV; // its VV relation
  vector<double> coords;   // vertex coordinates
  vector<uint> tris;       // triangle indices
  vector<double> res;

  vector<double> ground_truth;
  vector<double> estimated_distances;
  double smape;

  //-------- GUI state --------
  // View
  bool SHOW_MESH, SHOW_WIREFRAME;

  // Wireframe
  float wireframe_width; // Wireframe width
  float wireframe_alpha; // Wireframe transparency

  // Shading
  enum MeshRenderMode { RENDER_POINTS, RENDER_FLAT, RENDER_SMOOTH } render_mode;

  // Vector field
  DrawableVectorField vec_field;
  bool show_vecfield;
  float vecfield_size;
  Color vec_color;

  // Colors
  Color vert_color;
  Color poly_color;

  // Width of the Inputs
  float width;

  // SSGD Method
  enum SSGDMethod {
    VTP,
    TRETTNER,
    FAST_MARCHING,
    HEAT,
    GEOTANGLE,
    EDGE,
    EXTENDED,
    LANTHIER
  } ssgd_method;

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
  // // vector<int> vertices_idx = {1, 651, 1301, 1951, 2601};
  // vector<int> sources;

  // TODO: DA TOGLIERE
  // vector<int> sources = {93};
  int inputVertexIndex; // to store user input for vertex index
  vector<int> sources;

  // Trettner
  string mesh_path;
  HalfEdge mesh;

  VTPSolver vtp_solver;
  TrettnerSolver trettner_solver;
  FastMarchingSolver fast_mar_solver;
  HeatSolver heat_solver;
  GeotangleSolver geotangle_solver;
  EdgeSolver edge_solver;
  ExtendedSolver extended_solver;
  LanthierSolver lanthier_solver;

  // Timer
  std::chrono::steady_clock::time_point tic;
  std::chrono::steady_clock::time_point toc;
  double vtp_load, vtp_preprocess, vtp_query;
  double trettner_load, trettner_preprocess, trettner_query;
  double fast_mar_load, fast_mar_preprocess, fast_mar_query;
  double heat_load, heat_preprocess, heat_query;
  double geotangle_load, geotangle_preprocess, geotangle_query;
  double edge_load, edge_preprocess, edge_query;
  double extended_load, extended_preprocess, extended_query;
  double lanthier_load, lanthier_preprocess, lanthier_query;

  double true_FMM_query_time;

  State() {
    MESH_IS_LOADED = false;
    // View
    SHOW_WIREFRAME = false;
    SHOW_MESH = true;

    wireframe_width = 1.0;
    wireframe_alpha = 1.0;

    // Vector field
    show_vecfield = false;
    vecfield_size = 0.9f;
    vec_color = Color::RED();

    vert_color = Color::WHITE();
    poly_color = Color::WHITE();

    width = 100.0f;

    render_mode = State::RENDER_SMOOTH;
    ssgd_method = State::VTP;

    // Timer
    vtp_load, vtp_preprocess, vtp_query = 0.0;
    trettner_load, trettner_preprocess, trettner_query = 0.0;
    fast_mar_load, fast_mar_preprocess, fast_mar_query = 0.0;
    heat_load, heat_preprocess, heat_query = 0.0;
    geotangle_load, geotangle_preprocess, geotangle_query = 0.0;
    edge_load, edge_preprocess, edge_query = 0.0;
    extended_load, extended_preprocess, extended_query = 0.0;
    lanthier_load, lanthier_preprocess, lanthier_query = 0.0;

    true_FMM_query_time = 0.0;

    res = vector<double>();

    ground_truth = vector<double>();
    estimated_distances = vector<double>();
    smape = 0.0;

    mesh_path = "";
    mesh = HalfEdge();

    k = 3;
    prev_k = 3;

    n_steiner = 3;
    prev_n_steiner = 3;

    heat_time = 1.0;
    heat_time_prev = 1.0;
  }
};

void fillTimeTable(State &gs, const string &method_name, double load_time,
                   double preprocess_time, double query_time) {
  if (method_name == "VTP") {
    gs.vtp_load = load_time;
    // gs.vtp_preprocess = preprocess_time;
    gs.vtp_preprocess = 0.0;
    gs.vtp_query = query_time;
  } else if (method_name == "Trettner") {
    gs.trettner_load = load_time;
    // gs.trettner_preprocess = preprocess_time;
    gs.trettner_preprocess = 0.0;
    gs.trettner_query = query_time;
  } else if (method_name == "Fast Marching") {
    gs.fast_mar_load = load_time;
    // gs.fast_mar_preprocess = preprocess_time;
    gs.fast_mar_preprocess = 0.0;
    gs.fast_mar_query = query_time;
  } else if (method_name == "Heat") {
    gs.heat_load = load_time;
    gs.heat_preprocess = preprocess_time;
    gs.heat_query = query_time;
  } else if (method_name == "Geotangle") {
    gs.geotangle_load = load_time;
    gs.geotangle_preprocess = preprocess_time;
    gs.geotangle_query = query_time;
  } else if (method_name == "Edge") {
    gs.edge_load = load_time;
    gs.edge_preprocess = preprocess_time;
    gs.edge_query = query_time;
  } else if (method_name == "Extended") {
    gs.extended_load = load_time;
    gs.extended_preprocess = preprocess_time;
    gs.extended_query = query_time;
  } else if (method_name == "Lanthier") {
    gs.lanthier_load = load_time;
    gs.lanthier_preprocess = preprocess_time;
    gs.lanthier_query = query_time;
  }
}

void init(GeodesicMethod &m, State &gs, string name) {
  cout << "---------- Initializing method: " << name << " ----------" << endl;

  // load
  gs.tic = std::chrono::steady_clock::now();
  m.load(&gs.m);
  gs.toc = std::chrono::steady_clock::now();
  double load =
      chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic).count();
  cout << "Load time: " << load << " milliseconds" << endl;

  // preprocess
  // gs.tic = std::chrono::steady_clock::now();
  auto start = std::chrono::steady_clock::now();
  m.preprocess();
  // gs.toc = std::chrono::steady_clock::now();
  auto end = std::chrono::steady_clock::now();
  auto preprocess =
      std::chrono::duration_cast<std::chrono::duration<double>>(end - start)
          .count();
  cout << "Preprocess time: " << preprocess << " seconds" << endl;

  // fill time table
  fillTimeTable(gs, name, load, preprocess, 0.0);
}

void init_methods(State &gs, atomic<float> &progress) {
  init(gs.vtp_solver, gs, "VTP");
  progress.store(1.0f / 8.0f);

  // init(gs.trettner_solver, gs, "Trettner");
  // progress.store(2.0f / 8.0f);

  // init(gs.fast_mar_solver, gs, "Fast Marching");
  // progress.store(3.0f / 8.0f);

  // init(gs.heat_solver, gs, "Heat");
  // progress.store(4.0f / 8.0f);

  // init(gs.geotangle_solver, gs, "Geotangle");
  // progress.store(5.0f / 8.0f);

  // init(gs.edge_solver, gs, "Edge");
  // progress.store(6.0f / 8.0f);

  // init(gs.lanthier_solver, gs, "Lanthier");
  // progress.store(7.0f / 8.0f);

  // init(gs.extended_solver, gs, "Extended");
  // progress.store(8.0f / 8.0f);
}

// SMAPE calculation
// double calculate_smape(const vector<double>& gt, const vector<double>& est) {
//     cout << "GT size: " << gt.size() << ", EST size: " << est.size() << endl;
//     if (gt.empty() || est.empty()) {
//         cerr << "Ground truth or estimated distances are empty." << endl;
//         return 0.0;
//     }

//     double smape = 0.0;
//     int count = 0;

//     // for (size_t i = 0; i < gt.size() && i < est.size(); ++i) {
//     for (size_t i = 0; i < est.size(); ++i) {
//         double denom = std::abs(gt[i]) + std::abs(est[i]);
//         if (denom != 0) {
//             smape += std::abs(gt[i] - est[i]) / denom;
//             ++count;
//         }
//     }

//     if (count > 0) {
//         smape = (smape / count) * 100.0;  // Convert to percentage
//     }
//     return smape;
// }

// Updated SMAPE calculation function
std::pair<double, std::vector<double>>
calculate_smape(const std::vector<double> &gt, const std::vector<double> &est) {
  std::cout << "GT size: " << gt.size() << ", EST size: " << est.size()
            << std::endl;
  if (gt.empty() || est.empty()) {
    std::cerr << "Ground truth or estimated distances are empty." << std::endl;
    return {0.0, std::vector<double>()};
  }

  double smape_percentage = 0.0;
  int count = 0;
  std::vector<double> smape_values(gt.size(),
                                   0.0); // To store individual SMAPE values

  // Calculate SMAPE for each vertex and overall SMAPE percentage
  for (size_t i = 0; i < est.size(); ++i) {
    double denom = std::abs(gt[i]) + std::abs(est[i]);
    if (denom != 0) {
      smape_values[i] = std::abs(gt[i] - est[i]) / denom;
      smape_percentage += smape_values[i];
      ++count;
    }
  }

  // Convert to overall SMAPE percentage
  if (count > 0) {
    smape_percentage = (smape_percentage / count) * 100.0;
  }

  return {smape_percentage, smape_values};
}

Color interpolate_color(double value) {
  // Ensure the value is clamped between 0 and 1
  value = std::max(0.0, std::min(1.0, value));

  // Directly assign green to zero error without any transformation
  if (value == 0) {
    return Color(0.0f, 1.0f, 0.0f, 1.0f); // Pure green for zero error
  }

  // Apply a transformation to enhance contrast for non-zero values
  value = pow(value, 0.5); // Using square root to enhance contrast

  // Define the colors using RGBA format
  Color green(0.0f, 1.0f, 0.0f, 1.0f);  // Green color for low error
  Color yellow(1.0f, 1.0f, 0.0f, 1.0f); // Yellow color for medium error
  Color red(1.0f, 0.0f, 0.0f, 1.0f);    // Red color for high error

  float r, g, b, a = 1.0f; // Set alpha to fully opaque
  if (value < 0.5) {
    // Interpolate between green and yellow
    float local_value = value * 2; // Scale to 0-1 range
    r = (1.0f - local_value) * green.rgba[0] + local_value * yellow.rgba[0];
    g = (1.0f - local_value) * green.rgba[1] + local_value * yellow.rgba[1];
    b = (1.0f - local_value) * green.rgba[2] + local_value * yellow.rgba[2];
  } else {
    // Interpolate between yellow and red
    float local_value = (value - 0.5f) * 2; // Adjust scale for upper half
    r = (1.0f - local_value) * yellow.rgba[0] + local_value * red.rgba[0];
    g = (1.0f - local_value) * yellow.rgba[1] + local_value * red.rgba[1];
    b = (1.0f - local_value) * yellow.rgba[2] + local_value * red.rgba[2];
  }

  return Color(r, g, b, a);
}

void visualize_smape_on_mesh(State &gs,
                             const std::vector<double> &smape_values) {
  // Create a scalar field from the SMAPE values
  gs.field = ScalarField(smape_values);

  // Normalize the scalar field between 0 and 1 for better visualization
  gs.field.normalize_in_01();

  // Assign colors based on the normalized values using the red-blue gradient
  for (uint vid = 0; vid < gs.m.num_verts(); ++vid) {
    double normalized_value = gs.field[vid]; // Normalized SMAPE value

    // Use the custom interpolation function to get the color
    Color color = interpolate_color(normalized_value);

    // Set the interpolated color to the vertex
    gs.m.vert_data(vid).color = color;
  }

  // Show the vertex colors on the mesh
  gs.m.show_vert_color();
  gs.m.updateGL();
}

//:::::::::::::::::::::::::::::::::::: I/O ::::::::::::::::::::::::::::::::::::
void Load_mesh(string filename, GLcanvas &gui, State &gs) {
  gs.m = DrawableTrimesh<>(filename.c_str());
  gs.nverts = gs.m.num_verts();
  gs.VV.resize(gs.nverts); // fill in Vertex-Vertex relation

  gs.coords = extract_coords(gs.m);
  gs.tris = extract_tris(gs.m);

  for (auto i = 0; i < gs.nverts; i++)
    gs.VV[i] = gs.m.vert_ordered_verts_link(i);

  gs.m.normalize_bbox(); // rescale mesh to fit [0,1]^3 box
  gs.m.center_bbox();

  gs.m.show_wireframe(gs.SHOW_WIREFRAME);
  gs.m.show_mesh(gs.SHOW_MESH);
  gs.m.updateGL();

  if (gs.MESH_IS_LOADED) {

    // Clear and reinitialize the vector field for the new mesh
    gs.vec_field = DrawableVectorField();
    gs.show_vecfield = false; // Reset the flag to not show the old vector field
    // Clear the sources for the new mesh
    gs.sources_heat.clear(); // Reset the sources for the new mesh
    gs.sources.clear();      // Reset the sources for the new mesh
    // Clear cache for Heat method
    gs.prefactored_matrices.heat_flow_cache = NULL; // Reset the heat flow cache

    gs.vtp_load, gs.vtp_preprocess, gs.vtp_query = 0.0;
    gs.trettner_load, gs.trettner_preprocess, gs.trettner_query = 0.0;
    gs.fast_mar_load, gs.fast_mar_preprocess, gs.fast_mar_query = 0.0;
    gs.heat_load, gs.heat_preprocess, gs.heat_query = 0.0;
    gs.geotangle_load, gs.geotangle_preprocess, gs.geotangle_query = 0.0;
    gs.edge_load, gs.edge_preprocess, gs.edge_query = 0.0;
    gs.extended_load, gs.extended_preprocess, gs.extended_query = 0.0;
    gs.lanthier_load, gs.lanthier_preprocess, gs.lanthier_query = 0.0;

    gs.true_FMM_query_time = 0.0;

    gs.res = vector<double>();

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
    string out_normalized_bb_mesh =
        new_folder + "/" + new_filename.substr(0, new_filename.size() - 4) +
        "_NORMALIZED.obj";
    // Save the normalized mesh
    gs.m.save(out_normalized_bb_mesh.c_str());
    // Give to Trettner the normalized mesh
    gs.trettner_solver = TrettnerSolver(out_normalized_bb_mesh);

    // gs.mesh_path = filename;
    // gs.trettner_solver = TrettnerSolver(gs.mesh_path);

    init_methods(gs, progress);
  }

  if (!gs.MESH_IS_LOADED) {
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
    string out_normalized_bb_mesh =
        new_folder + "/" + new_filename.substr(0, new_filename.size() - 4) +
        "_NORMALIZED.obj";
    // Save the normalized mesh
    gs.m.save(out_normalized_bb_mesh.c_str());
    // Give to Trettner the normalized mesh
    gs.trettner_solver = TrettnerSolver(out_normalized_bb_mesh);

    // gs.mesh_path = filename;
    // gs.trettner_solver = TrettnerSolver(gs.mesh_path);

    gui.push(&gs.m);
    gs.MESH_IS_LOADED = true;
  }
}

void Load_mesh(GLcanvas &gui, State &gs) {
  string filename = file_dialog_open();
  if (filename.size() != 0)
    gs.mesh_path = filename;

  gs.mesh_path = filename;
  Load_mesh(filename, gui, gs);
}

//::::::::::::::::::::::::::::::::::::::
//: GUI:::::::::::::::::::::::::::::::::::::::::::::::::
GLcanvas Init_GUI() {
  GLcanvas gui(1600, 800);
  gui.side_bar_width = 0.30;
  gui.show_side_bar = true;
  return gui;
}

void Setup_GUI_Callbacks(GLcanvas &gui, State &gs) {
  gui.callback_app_controls = [&]() {
    // New detached window
    bool show_new_window = true;

    if (show_new_window) {
      float sidebar_width = gui.side_bar_width * 1500;
      ImVec2 new_window_pos = ImVec2(sidebar_width + 650, 25);
      ImVec2 window_size = ImVec2(400, 200);

      ImGui::SetNextWindowPos(new_window_pos, ImGuiCond_FirstUseEver);
      ImGui::SetNextWindowSize(window_size, ImGuiCond_FirstUseEver);

      ImGui::Begin("SSGD Methods Timing Results", &show_new_window);

      // Define the table headers
      if (ImGui::BeginTable("TimingResults", 4)) {
        ImGui::TableSetupColumn("Method");
        ImGui::TableSetupColumn("Load (ms)");
        ImGui::TableSetupColumn("Preprocess (ms)");
        ImGui::TableSetupColumn("Query (ms)");
        ImGui::TableHeadersRow();

        // VTP Method
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::Text("VTP");
        ImGui::TableSetColumnIndex(1);
        ImGui::Text("%.4f", gs.vtp_load);
        ImGui::TableSetColumnIndex(2);
        ImGui::Text("%.4f", gs.vtp_preprocess);
        ImGui::TableSetColumnIndex(3);
        ImGui::Text("%.4f", gs.vtp_query);

        // Trettner Method
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::Text("Trettner");
        ImGui::TableSetColumnIndex(1);
        ImGui::Text("%.4f", gs.trettner_load);
        ImGui::TableSetColumnIndex(2);
        ImGui::Text("%.4f", gs.trettner_preprocess);
        ImGui::TableSetColumnIndex(3);
        ImGui::Text("%.4f", gs.trettner_query);

        // Fast Marching Method
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::Text("Fast Marching");
        ImGui::TableSetColumnIndex(1);
        ImGui::Text("%.4f", gs.fast_mar_load);
        ImGui::TableSetColumnIndex(2);
        ImGui::Text("%.4f", gs.fast_mar_preprocess);
        ImGui::TableSetColumnIndex(3);
        ImGui::Text("%.4f", gs.fast_mar_query);

        // Heat Method
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::Text("Heat");
        ImGui::TableSetColumnIndex(1);
        ImGui::Text("%.4f", gs.heat_load);
        ImGui::TableSetColumnIndex(2);
        ImGui::Text("%.4f", gs.heat_preprocess);
        ImGui::TableSetColumnIndex(3);
        ImGui::Text("%.4f", gs.heat_query);

        // GeoTangle Method
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::Text("GeoTangle");
        ImGui::TableSetColumnIndex(1);
        ImGui::Text("%.4f", gs.geotangle_load);
        ImGui::TableSetColumnIndex(2);
        ImGui::Text("%.4f", gs.geotangle_preprocess);
        ImGui::TableSetColumnIndex(3);
        ImGui::Text("%.4f", gs.geotangle_query);

        // Edge Method
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::Text("Edge");
        ImGui::TableSetColumnIndex(1);
        ImGui::Text("%.4f", gs.edge_load);
        ImGui::TableSetColumnIndex(2);
        ImGui::Text("%.4f", gs.edge_preprocess);
        ImGui::TableSetColumnIndex(3);
        ImGui::Text("%.4f", gs.edge_query);

        // Extended Method
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::Text("Extended");
        ImGui::TableSetColumnIndex(1);
        ImGui::Text("%.4f", gs.extended_load);
        ImGui::TableSetColumnIndex(2);
        ImGui::Text("%.4f", gs.extended_preprocess);
        ImGui::TableSetColumnIndex(3);
        ImGui::Text("%.4f", gs.extended_query);

        // Lanthier Method
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::Text("Lanthier");
        ImGui::TableSetColumnIndex(1);
        ImGui::Text("%.4f", gs.lanthier_load);
        ImGui::TableSetColumnIndex(2);
        ImGui::Text("%.4f", gs.lanthier_preprocess);
        ImGui::TableSetColumnIndex(3);
        ImGui::Text("%.4f", gs.lanthier_query);

        ImGui::EndTable();
      }

      ImGui::End();
    }

    // In your main function or GUI rendering loop
    if (progress < 1.0f) {
      ImVec2 window_size = ImVec2(250, 60);  // Width, Height in pixels
      ImVec2 window_pos = ImVec2(1225, 300); // X, Y position in pixels

      ImGui::SetNextWindowSize(window_size, ImGuiCond_Always);
      ImGui::SetNextWindowPos(window_pos, ImGuiCond_Always);

      ImGui::Begin("Graph Construction Progress");
      ImGui::ProgressBar(progress.load(), ImVec2(-1.0f, 0.0f),
                         progress > 1.0f ? "Done" : "Loading...");
      ImGui::End();
    }

    // Files
    ImGui::PushFont(lato_bold_title);
    ImGui::SeparatorText("Single Source Geodesic Distance Computation");
    ImGui::PopFont();

    ImGui::SetNextItemOpen(true, ImGuiCond_Once);
    ImGui::PushFont(lato_bold);
    if (ImGui::TreeNode("IO")) {
      ImGui::PopFont();
      if (ImGui::Button("Load mesh")) {
        if (gs.MESH_IS_LOADED) {
          ImGui::OpenPopup("Load mesh?");
        } else {
          Load_mesh(gui, gs);
        }
      }
      ImGui::SameLine();
      if (ImGui::SmallButton("Save")) {
        std::string filename = file_dialog_save();
        if (!filename.empty()) {
          gs.m.save(filename.c_str());
        }
      }

      ImGui::PushFont(lato_bold);
      ImGui::SeparatorText("Mesh Information");
      ImGui::PopFont();
      // Assuming Load_mesh successfully loads the mesh into gs.m
      int numVertices = gs.m.num_verts();
      int numFaces =
          gs.m.num_polys(); // or num_faces() depending on your mesh type

      ImGui::Text("Number of vertices: %d", numVertices);
      ImGui::Text("Number of faces: %d", numFaces);
      ImGui::Text("");

      // Modal popup for loading files
      ImVec2 center = ImGui::GetMainViewport()->GetCenter();
      ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
      if (ImGui::BeginPopupModal("Load mesh?", NULL,
                                 ImGuiWindowFlags_AlwaysAutoResize |
                                     ImGuiWindowFlags_NoCollapse)) {
        static bool dont_ask_me_next_time = false;

        if (dont_ask_me_next_time) {
          Load_mesh(gui, gs);
          ImGui::CloseCurrentPopup();
        }
        ImGui::Text("All data structures will be reset - Load anyway?\n\n");
        ImGui::Separator();
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(0, 0));
        ImGui::Checkbox("Don't ask me next time", &dont_ask_me_next_time);
        ImGui::PopStyleVar();
        if (ImGui::Button("OK", ImVec2(120, 0))) {
          Load_mesh(gui, gs);
          ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel", ImVec2(120, 0))) {
          ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
      }
      ImGui::TreePop();
    } else {
      ImGui::PopFont();
    }

    ImGui::Text("Set Source Vertex:");
    ImGui::InputInt("Vertex Index", &gs.inputVertexIndex);
    if (ImGui::Button("Set")) {
      if (gs.inputVertexIndex >= 0 && gs.inputVertexIndex < gs.m.num_verts()) {
        gs.sources.clear(); // Clear existing sources if you only want one
                            // source at a time
        gs.sources.push_back(gs.inputVertexIndex);
        std::cout << "Source vertex set to: " << gs.inputVertexIndex
                  << std::endl;
      } else {
        std::cerr << "Invalid vertex index!" << std::endl;
      }
    }
    ImGui::Separator();

    // Wireframe settings
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    ImGui::PushFont(lato_bold);
    if (ImGui::TreeNode("Wireframe")) {
      ImGui::PopFont();
      if (ImGui::Checkbox("Show", &gs.SHOW_WIREFRAME))
        gs.m.show_wireframe(gs.SHOW_WIREFRAME);
      if (ImGui::SliderFloat("Width", &gs.wireframe_width, 1.f, 10.f))
        gs.m.show_wireframe_width(gs.wireframe_width);
      if (ImGui::SliderFloat("Transparency", &gs.wireframe_alpha, 0.f, 1.f))
        gs.m.show_wireframe_transparency(gs.wireframe_alpha);
      ImGui::Text("");
      ImGui::TreePop();
    } else {
      ImGui::PopFont();
    }

    // Mesh shading
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    ImGui::PushFont(lato_bold);
    if (ImGui::TreeNode("Shading")) {
      ImGui::PopFont();
      // Point rendering mode
      if (ImGui::RadioButton("Point ",
                             gs.render_mode == State::RENDER_POINTS)) {
        gs.render_mode = State::RENDER_POINTS;
        gs.m.show_mesh_points();
      }
      // Flat shading mode
      if (ImGui::RadioButton("Flat  ", gs.render_mode == State::RENDER_FLAT)) {
        gs.render_mode = State::RENDER_FLAT;
        gs.m.show_mesh_flat();
      }
      // Smooth shading mode
      if (ImGui::RadioButton("Smooth",
                             gs.render_mode == State::RENDER_SMOOTH)) {
        gs.render_mode = State::RENDER_SMOOTH;
        gs.m.show_mesh_smooth(); // Assuming you have a method to display smooth
      }
      ImGui::Text("");
      ImGui::TreePop();
    } else {
      ImGui::PopFont();
    }

    // Mesh colors
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    ImGui::PushFont(lato_bold);
    if (ImGui::TreeNode("Colors")) {
      ImGui::PopFont();
      if (ImGui::BeginTable("Color by:", 2)) {
        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        if (ImGui::RadioButton("Vert",
                               gs.m.drawlist.draw_mode & DRAW_TRI_VERTCOLOR)) {
          gs.m.show_vert_color();
        }
        ImGui::TableNextColumn();
        if (ImGui::ColorEdit4("Vertex Color", gs.vert_color.rgba)) {
          gs.m.vert_set_color(gs.vert_color);
          gs.m.show_vert_color();
          gs.m.updateGL();
        }

        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        if (ImGui::RadioButton("Poly",
                               gs.m.drawlist.draw_mode & DRAW_TRI_FACECOLOR)) {
          gs.m.show_poly_color();
        }
        ImGui::TableNextColumn();
        if (ImGui::ColorEdit4("Polygon Color", gs.poly_color.rgba)) {
          gs.m.poly_set_color(gs.poly_color);
          gs.m.show_poly_color();
          gs.m.updateGL();
        }

        ImGui::EndTable();
      }
      ImGui::Text("");
      ImGui::TreePop();
    } else {
      ImGui::PopFont();
    }

    // Vector field visualization
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    ImGui::PushFont(lato_bold);
    if (ImGui::TreeNode("Vector Field")) {
      ImGui::PopFont();
      if (ImGui::Checkbox("Show Vector Field", &gs.show_vecfield)) {
        if (gs.show_vecfield) {
          if (gs.vec_field.size() ==
              0) { // Initialize the vector field if not already
            gs.vec_field = DrawableVectorField(gs.m);
            ScalarField f(gs.m.serialize_uvw(3)); // Adjust U_param if needed
            gs.vec_field = gradient_matrix(gs.m) * f;
            gs.vec_field.normalize();
            // print the values inside vec_field
            // Print the number of elements in the vector field
            gs.vec_field.set_arrow_size(float(gs.m.edge_avg_length()) *
                                        gs.vecfield_size);
            gs.vec_field.set_arrow_color(gs.vec_color);
            gs.m.updateGL();
          }
          gui.push(&gs.vec_field, false);
        } else {
          gui.pop(&gs.vec_field);
          gs.m.updateGL();
        }
      }

      // Vector field size control
      if (ImGui::SliderFloat("Size", &gs.vecfield_size, 0.1f, 5.f)) {
        gs.vec_field.set_arrow_size(float(gs.m.edge_avg_length()) *
                                    gs.vecfield_size);
      }

      // Vector field color control
      if (ImGui::ColorEdit4("Color##vec", gs.vec_color.rgba)) {
        gs.vec_field.set_arrow_color(gs.vec_color);
      }
      ImGui::Text("");
      ImGui::TreePop();
    } else {
      ImGui::PopFont();
    }

    // SSGD Method
    ImGui::SetNextItemOpen(true, ImGuiCond_Once);
    ImGui::PushFont(lato_bold);
    if (ImGui::TreeNode("SSGD Method")) {
      ImGui::PopFont();

      // Exact Polyhedral Methods
      ImGui::PushFont(lato_bold);
      ImGui::SeparatorText("Exact Polyhedral Methods");
      ImGui::PopFont();
      if (ImGui::RadioButton("VTP ", gs.ssgd_method == State::VTP)) {
        gs.ssgd_method = State::VTP;
      }
      if (ImGui::RadioButton("Trettner ", gs.ssgd_method == State::TRETTNER)) {
        gs.ssgd_method = State::TRETTNER;
      }

      // PDE-Based Methods
      ImGui::PushFont(lato_bold);
      ImGui::SeparatorText("PDE-Based Methods");
      ImGui::PopFont();

      if (ImGui::RadioButton("Fast Marching  ",
                             gs.ssgd_method == State::FAST_MARCHING)) {
        gs.ssgd_method = State::FAST_MARCHING;
      }

      // Define two columns
      ImGui::Columns(2, nullptr, false);
      // First column for the radio button
      if (ImGui::RadioButton("Heat  ", gs.ssgd_method == State::HEAT)) {
        gs.ssgd_method = State::HEAT;
      }
      ImGui::NextColumn();
      // Second column for the InputDouble
      ImGui::SetNextItemWidth(gs.width);
      ImGui::InputDouble("t", &gs.heat_time, 0.1, 1.0, "%.3f");
      ImGui::Columns(1);

      // Graph-Based Methods
      ImGui::PushFont(lato_bold);
      ImGui::SeparatorText("Graph-Based Methods");
      ImGui::PopFont();
      if (ImGui::RadioButton("GeoTangle  ",
                             gs.ssgd_method == State::GEOTANGLE)) {
        gs.ssgd_method = State::GEOTANGLE;
      }
      if (ImGui::RadioButton("Edge  ", gs.ssgd_method == State::EDGE)) {
        gs.ssgd_method = State::EDGE;
      }

      ImGui::Columns(2, nullptr, false);
      // First column for the radio button
      if (ImGui::RadioButton("Extended  ", gs.ssgd_method == State::EXTENDED)) {
        gs.ssgd_method = State::EXTENDED;
      }
      ImGui::NextColumn();
      // Second column for the InputInt
      ImGui::SetNextItemWidth(gs.width);
      ImGui::InputInt("k", &gs.k, 1,
                      10); // Add an input integer for parameter k
      ImGui::Columns(1);
      gs.k = std::max(gs.k, 1);

      ImGui::Columns(2, nullptr, false);
      if (ImGui::RadioButton("Lanthier  ", gs.ssgd_method == State::LANTHIER)) {
        gs.ssgd_method = State::LANTHIER;
      }
      ImGui::NextColumn();
      // Second column for the InputInt
      ImGui::SetNextItemWidth(gs.width);
      ImGui::InputInt("Steiner points", &gs.n_steiner, 1, 10);
      ImGui::Columns(1);
      gs.n_steiner = std::max(gs.n_steiner, 1);
      // cout << "Number of Steiner points: " << gs.n_steiner << endl;

      ImGui::TreePop();
    } else {
      ImGui::PopFont();
    }

    ImGui::Text("");
    ImGui::Text("Pick a source vertex: shift + left click");

    // Display the label "Sources"
    ImGui::Text("Sources:");
    // Iterate over the gs.sources to display each source vertex ID
    for (uint i = 0; i < gs.sources_heat.size(); ++i) {
      ImGui::Text("Vertex ID: %u", gs.sources_heat[i]);
    }
    ImGui::Text("");

    // Button for Compute SSGD
    if (ImGui::Button("Compute SSGD")) {
      // Based on the selected SSGD method, perform different actions
      if (gs.sources_heat.empty() && gs.sources.empty()) {
        // Open a warning popup if no source is selected
        ImGui::OpenPopup("Warning");
      } else {

        switch (gs.ssgd_method) {

        case State::VTP: {
          // load() and preprocess() in init_methods()
          gs.tic = std::chrono::steady_clock::now();
          gs.vtp_solver.query(gs.sources[0], gs.res);
          gs.toc = std::chrono::steady_clock::now();
          gs.vtp_query = chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic).count();
          // for (int i = 0; i < gs.res.size(); i++) {
          //   cout << i << "," << gs.res[i] << endl;
          // }
          // cout << "---------------" << endl;
          // cout << "VTP: " << endl;
          // cout << "---------------" << endl;
          // cout << "Source: " << gs.sources[0] << endl;
          // for (int i = 0; i < gs.res.size(); i++) {
          //   // cout << i << "," << gs.res[i] << endl;
          //   cout << gs.res[i] << "," << endl;
          // }
          double seconds = std::chrono::duration_cast<std::chrono::duration<double>>(gs.toc - gs.tic).count();
          cout << "QUERY TIME: " << seconds << " seconds" << std::endl;

          // I want to write gs.res in a .txt file one value at row
          // Create the file in this directory 
  
          

          fillTimeTable(gs, "VTP", gs.vtp_load, gs.vtp_preprocess, gs.vtp_query);
          break;
        }

        case State::TRETTNER: {
          gs.tic = std::chrono::steady_clock::now();
          gs.trettner_solver.query(gs.sources[0], gs.res);
          gs.toc = std::chrono::steady_clock::now();
          gs.trettner_query =
              chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic)
                  .count();
          // cout << "---------------" << endl;
          // cout << "TRETTNER: " << endl;
          // cout << "---------------" << endl;
          // for (int i = 0; i < gs.res.size(); i++) {
          //   cout << i << "," << gs.res[i] << endl;
          // }
          double seconds =
              std::chrono::duration_cast<std::chrono::duration<double>>(gs.toc -
                                                                        gs.tic)
                  .count();
          cout << "QUERY TIME: " << seconds << " seconds" << std::endl;

          fillTimeTable(gs, "Trettner", gs.trettner_load,
                        gs.trettner_preprocess, gs.trettner_query);
          break;
        }

        case State::FAST_MARCHING: {
          gs.tic = std::chrono::steady_clock::now();
          gs.fast_mar_solver.query(gs.sources[0], gs.res);
          gs.toc = std::chrono::steady_clock::now();
          gs.fast_mar_query =
              chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic)
                  .count();
          // cout << "---------------" << endl;
          // cout << "FAST MAR: " << endl;
          // cout << "---------------" << endl;
          // for (int i = 0; i < gs.res.size(); i++) {
          //   cout << i << "," << gs.res[i] << endl;
          // }

          // take query time
          gs.true_FMM_query_time = gs.fast_mar_solver.get_time();
          cout << "Query function time: " << gs.fast_mar_query << endl;
          cout << "True FMM query time: " << gs.true_FMM_query_time << endl;

          fillTimeTable(gs, "Fast Marching", gs.fast_mar_load,
                        gs.fast_mar_preprocess, gs.true_FMM_query_time);
          break;
        }

        case State::HEAT: {
          if (gs.heat_time != gs.heat_time_prev) {
            cout << "Heat time has changed to: " << gs.heat_time << endl;
            gs.heat_solver.set_t(gs.heat_time); // set_t recall the preprocess()
            gs.heat_time_prev = gs.heat_time;
          }

          double time_scalar = 1;

          cout << "Time scalar: " << time_scalar << endl;
          gs.heat_solver.set_t(time_scalar);

          gs.tic = std::chrono::steady_clock::now();
          gs.heat_solver.query(gs.sources[0], gs.res);
          gs.toc = std::chrono::steady_clock::now();
          gs.heat_query =
              chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic)
                  .count();
          // cout << "---------------" << endl;
          // cout << "HEAT: " << endl;
          // cout << "---------------" << endl;
          // for (int i = 0; i < gs.res.size(); i++) {
          //   cout << i << "," << gs.res[i] << endl;
          // }
          double seconds =
              std::chrono::duration_cast<std::chrono::duration<double>>(gs.toc -
                                                                        gs.tic)
                  .count();
          cout << "QUERY TIME: " << seconds << " seconds" << std::endl;
          ;

          fillTimeTable(gs, "Heat", gs.heat_load, gs.heat_preprocess,
                        gs.heat_query);
          break;
        }

        case State::GEOTANGLE: {
          gs.tic = std::chrono::steady_clock::now();
          gs.geotangle_solver.query(gs.sources[0], gs.res);
          gs.toc = std::chrono::steady_clock::now();
          gs.geotangle_query =
              chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic)
                  .count();

          // cout << "---------------" << endl;
          // cout << "GEOTANGLE: " << endl;
          // cout << "---------------" << endl;
          // for (int i = 0; i < gs.res.size(); i++) {
          //   cout << i << "," << gs.res[i] << endl;
          // }
          double seconds =
              std::chrono::duration_cast<std::chrono::duration<double>>(gs.toc -
                                                                        gs.tic)
                  .count();
          cout << "QUERY TIME: " << seconds << " seconds" << std::endl;

          fillTimeTable(gs, "Geotangle", gs.geotangle_load,
                        gs.geotangle_preprocess, gs.geotangle_query);
          break;
        }

        case State::EDGE: {
          gs.tic = std::chrono::steady_clock::now();
          gs.edge_solver.query(gs.sources[0], gs.res);

          // cout << "---------------" << endl;
          // cout << "EDGE: " << endl;
          // cout << "---------------" << endl;
          // for (int i = 0; i < gs.res.size(); i++) {
          //   cout << gs.res[i] << endl;
          // }
          gs.toc = std::chrono::steady_clock::now();
          gs.edge_query =
              chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic)
                  .count();

          double seconds =
              std::chrono::duration_cast<std::chrono::duration<double>>(gs.toc -
                                                                        gs.tic)
                  .count();
          cout << "QUERY TIME: " << seconds << " seconds" << std::endl;

          fillTimeTable(gs, "Edge", gs.edge_load, gs.edge_preprocess,
                        gs.edge_query);
          break;
        }

        case State::EXTENDED: {
          if (gs.k != gs.prev_k) {
            cout << "k has changed to: " << gs.k << endl;

            gs.tic = std::chrono::steady_clock::now();
            gs.extended_solver.set_k(gs.k); // This is the blocking call
            gs.toc = std::chrono::steady_clock::now();

            // Calculate the time taken for preprocessing
            gs.extended_preprocess =
                chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic)
                    .count();
            fillTimeTable(gs, "Extended", gs.extended_load,
                          gs.extended_preprocess, gs.extended_query);

            gs.prev_k = gs.k;
          }

          gs.tic = std::chrono::steady_clock::now();
          gs.extended_solver.query(gs.sources[0], gs.res);
          // cout << "---------------" << endl;
          // cout << "EXTENDED: " << endl;
          // cout << "---------------" << endl;
          // for (int i = 0; i < gs.res.size(); i++) {
          //   cout << gs.res[i] << endl;
          // }
          gs.toc = std::chrono::steady_clock::now();
          gs.extended_query =
              chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic)
                  .count();

          double seconds =
              std::chrono::duration_cast<std::chrono::duration<double>>(gs.toc -
                                                                        gs.tic)
                  .count();
          cout << "QUERY TIME: " << seconds << " seconds" << std::endl;

          fillTimeTable(gs, "Extended", gs.extended_load,
                        gs.extended_preprocess, gs.extended_query);
          break;
        }

        case State::LANTHIER: {
          if (gs.n_steiner != gs.prev_n_steiner) {
            cout << "Number of Steiner points has changed to: " << gs.n_steiner
                 << endl;

            gs.tic = std::chrono::steady_clock::now();
            gs.lanthier_solver.set_n_steiner(gs.n_steiner);
            gs.toc = std::chrono::steady_clock::now();

            // Calculate the time taken for preprocessing
            gs.lanthier_preprocess =
                chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic)
                    .count();
            fillTimeTable(gs, "Lanthier", gs.lanthier_load,
                          gs.lanthier_preprocess, gs.lanthier_query);

            gs.prev_n_steiner = gs.n_steiner;
          }

          gs.tic = std::chrono::steady_clock::now();
          cout << "Source vertex: " << gs.sources[0] << endl;
          gs.lanthier_solver.query(gs.sources[0], gs.res);
          // cout << "---------------" << endl;
          // cout << "LANTHIER: " << endl;
          // cout << "---------------" << endl;
          // for (int i = 0; i < gs.res.size(); i++) {
          //   cout << i << "," << gs.res[i] << endl;
          // }
          gs.toc = std::chrono::steady_clock::now();
          gs.lanthier_query =
              chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic)
                  .count();

          double seconds =
              std::chrono::duration_cast<std::chrono::duration<double>>(gs.toc -
                                                                        gs.tic)
                  .count();
          cout << "QUERY TIME: " << seconds << " seconds" << std::endl;

          fillTimeTable(gs, "Lanthier", gs.lanthier_load,
                        gs.lanthier_preprocess, gs.lanthier_query);
          break;
        }

        default:
          cout << "No SSGD method selected" << endl;
          break;
        }
        gs.field = ScalarField(gs.res);
        gs.field.normalize_in_01();
        gs.field.copy_to_mesh(gs.m);
        gs.m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);
      }
    }

    static float displayed_smape = 0.0f;
    // Button for Compute SMAPE error
    if (ImGui::Button("Compute SMAPE")) {
      // Based on the selected SSGD method, perform different actions
      if (gs.sources_heat.empty() && gs.sources.empty()) {
        // Open a warning popup if no source is selected
        ImGui::OpenPopup("Warning");
      } else {

        switch (gs.ssgd_method) {

        case State::VTP: {
          gs.vtp_solver.query(gs.sources[0], gs.res);
          gs.ground_truth = gs.res;
          gs.vtp_solver.query(gs.sources[0], gs.res);
          // gs.smape = calculate_smape(gs.ground_truth, gs.res);

          auto [smape_percentage, smape_values] =
              calculate_smape(gs.ground_truth, gs.res);
          gs.smape = smape_percentage;
          visualize_smape_on_mesh(gs, smape_values);

          cout << "SMAPE ERROR for VTP: " << gs.smape << endl;

          gs.ground_truth.clear();
          // gs.smape = 0.0;

          break;
        }

        case State::TRETTNER: {
          gs.res.clear();
          gs.vtp_solver.query(gs.sources[0], gs.res);
          gs.ground_truth = gs.res;
          gs.res.clear();
          gs.trettner_solver.load(&gs.m);
          gs.trettner_solver.preprocess();
          gs.trettner_solver.query(gs.sources[0], gs.res);
          // gs.smape = calculate_smape(gs.ground_truth, gs.res);

          auto [smape_percentage, smape_values] =
              calculate_smape(gs.ground_truth, gs.res);
          gs.smape = smape_percentage;
          visualize_smape_on_mesh(gs, smape_values);

          cout << "SMAPE ERROR for Trettner: " << gs.smape << endl;

          gs.ground_truth.clear();
          // gs.smape = 0.0;

          break;
        }

        case State::FAST_MARCHING: {
          gs.vtp_solver.query(gs.sources[0], gs.res);
          gs.ground_truth = gs.res;
          gs.fast_mar_solver.query(gs.sources[0], gs.res);
          // gs.smape = calculate_smape(gs.ground_truth, gs.res);

          auto [smape_percentage, smape_values] =
              calculate_smape(gs.ground_truth, gs.res);
          gs.smape = smape_percentage;
          visualize_smape_on_mesh(gs, smape_values);

          cout << "SMAPE ERROR for Fast Marching: " << gs.smape << endl;

          gs.ground_truth.clear();
          // gs.smape = 0.0;

          break;
        }

        case State::HEAT: {
          if (gs.heat_time != gs.heat_time_prev) {
            cout << "Heat time has changed to: " << gs.heat_time << endl;
            gs.heat_solver.set_t(gs.heat_time); // set_t recall the preprocess()
            gs.heat_time_prev = gs.heat_time;
          }

          double time_scalar = 1;

          cout << "Time scalar: " << time_scalar << endl;
          gs.heat_solver.set_t(time_scalar);

          gs.vtp_solver.query(gs.sources[0], gs.res);
          gs.ground_truth = gs.res;
          gs.heat_solver.query(gs.sources[0], gs.res);
          // gs.smape = calculate_smape(gs.ground_truth, gs.res);

          auto [smape_percentage, smape_values] =
              calculate_smape(gs.ground_truth, gs.res);
          gs.smape = smape_percentage;
          visualize_smape_on_mesh(gs, smape_values);

          cout << "SMAPE ERROR for Heat: " << gs.smape << endl;

          gs.ground_truth.clear();
          // gs.smape = 0.0;

          break;
        }

        case State::GEOTANGLE: {
          gs.vtp_solver.query(gs.sources[0], gs.res);
          gs.ground_truth = gs.res;
          gs.geotangle_solver.query(gs.sources[0], gs.res);
          // gs.smape = calculate_smape(gs.ground_truth, gs.res);

          auto [smape_percentage, smape_values] =
              calculate_smape(gs.ground_truth, gs.res);
          gs.smape = smape_percentage;
          visualize_smape_on_mesh(gs, smape_values);

          cout << "SMAPE ERROR for Geotangle: " << gs.smape << endl;

          gs.ground_truth.clear();
          // gs.smape = 0.0;

          break;
        }

        case State::EDGE: {
          gs.vtp_solver.query(gs.sources[0], gs.res);
          gs.ground_truth = gs.res;
          gs.edge_solver.query(gs.sources[0], gs.res);
          // gs.smape = calculate_smape(gs.ground_truth, gs.res);

          auto [smape_percentage, smape_values] =
              calculate_smape(gs.ground_truth, gs.res);
          gs.smape = smape_percentage;
          visualize_smape_on_mesh(gs, smape_values);

          cout << "SMAPE ERROR for Edge: " << gs.smape << endl;

          gs.ground_truth.clear();
          // gs.smape = 0.0;

          break;
        }

        case State::EXTENDED: {
          if (gs.k != gs.prev_k) {
            cout << "k has changed to: " << gs.k << endl;

            gs.tic = std::chrono::steady_clock::now();
            gs.extended_solver.set_k(gs.k); // This is the blocking call
            gs.toc = std::chrono::steady_clock::now();

            // Calculate the time taken for preprocessing
            gs.extended_preprocess =
                chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic)
                    .count();
            fillTimeTable(gs, "Extended", gs.extended_load,
                          gs.extended_preprocess, gs.extended_query);

            gs.prev_k = gs.k;
          }

          gs.vtp_solver.query(gs.sources[0], gs.res);
          gs.ground_truth = gs.res;
          gs.extended_solver.query(gs.sources[0], gs.res);
          // gs.smape = calculate_smape(gs.ground_truth, gs.res);

          auto [smape_percentage, smape_values] =
              calculate_smape(gs.ground_truth, gs.res);
          gs.smape = smape_percentage;
          visualize_smape_on_mesh(gs, smape_values);

          cout << "SMAPE ERROR for Extended: " << gs.smape << endl;

          gs.ground_truth.clear();
          // gs.smape = 0.0;

          break;
        }

        case State::LANTHIER: {
          if (gs.n_steiner != gs.prev_n_steiner) {
            cout << "Number of Steiner points has changed to: " << gs.n_steiner
                 << endl;

            gs.tic = std::chrono::steady_clock::now();
            gs.lanthier_solver.set_n_steiner(gs.n_steiner);
            gs.toc = std::chrono::steady_clock::now();

            // Calculate the time taken for preprocessing
            gs.lanthier_preprocess =
                chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic)
                    .count();
            fillTimeTable(gs, "Lanthier", gs.lanthier_load,
                          gs.lanthier_preprocess, gs.lanthier_query);

            gs.prev_n_steiner = gs.n_steiner;
          }

          gs.vtp_solver.query(gs.sources[0], gs.res);
          gs.ground_truth = gs.res;
          gs.lanthier_solver.query(gs.sources[0], gs.res);
          // gs.smape = calculate_smape(gs.ground_truth, gs.res);

          auto [smape_percentage, smape_values] =
              calculate_smape(gs.ground_truth, gs.res);
          gs.smape = smape_percentage;
          visualize_smape_on_mesh(gs, smape_values);

          cout << "SMAPE ERROR for Lanthier: " << gs.smape << endl;

          gs.ground_truth.clear();
          // gs.smape = 0.0;

          break;
        }

        default:
          cout << "No SSGD method selected" << endl;
          break;
        }
        // QUI VOLGIO VISUALIZZARE L'ERRORE COME UNA HEATMAP SULLA MESH PERCHE'
        // L'ERRORE E' UNO SCALAR gs.field = ScalarField(gs.res);
        // gs.field.normalize_in_01();
        // gs.field.copy_to_mesh(gs.m);
        // gs.m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);
      }
    }
    // Display SMAPE error next to the button
    ImGui::SameLine();
    ImGui::SetCursorPosX(ImGui::GetCursorPosX() + 20.0f);

    ImGui::PushFont(lato_bold);
    ImGui::Text("SMAPE Error:");
    ImGui::PopFont();

    ImGui::SameLine();
    ImGui::Text("%.6f%%", gs.smape);

    // Popup modal for warning if no source is selected
    if (ImGui::BeginPopupModal("Warning", NULL,
                               ImGuiWindowFlags_AlwaysAutoResize)) {
      ImGui::Text("Select a source vertex before computing SSGD");
      ImGui::Text("shift + left click");
      if (ImGui::Button("OK")) {
        ImGui::CloseCurrentPopup();
      }
      ImGui::EndPopup();
    }

    // Reset button
    if (ImGui::SmallButton("Reset")) {
      // Reset sources
      gs.sources_heat.clear();
      gs.sources.clear();

      // TODO: DA TOGLIEREEEEE
      // gs.sources = {93};

      gs.vtp_query = 0.0;
      gs.trettner_query = 0.0;
      gs.heat_query = 0.0;
      gs.geotangle_query = 0.0;
      gs.edge_query = 0.0;
      gs.extended_query = 0.0;
      gs.lanthier_query = 0.0;

      gs.res = vector<double>();

      // Reset the smape error
      gs.ground_truth = vector<double>();
      gs.estimated_distances = vector<double>();
      gs.smape = 0.0;

      // Reset the scalar field
      for (uint vid = 0; vid < gs.m.num_verts(); ++vid) {
        gs.m.vert_data(vid).color =
            Color::WHITE(); // Replace `original_color` with the actual color
      }
      gs.m.show_vert_color();
    }
  };
}

// Compute the geodesic distances from a single source vertex
void Setup_Mouse_Callback(GLcanvas &gui, State &gs) {
  gui.callback_mouse_left_click = [&](int modifiers) -> bool {
    if (modifiers & GLFW_MOD_SHIFT) {
      vec3d p;
      vec2d click = gui.cursor_pos();
      if (gui.unproject(click, p)) {
        // ----- HEAT SOURCES -----
        uint vid = gs.m.pick_vert(p);
        gs.sources_heat.push_back(vid);
        cout << "Source vertex: " << vid << endl;

        // ------ SOURCES ------
        int selected_vid = gs.m.pick_vert(p);
        gs.sources.push_back(selected_vid);

        gs.m.vert_data(vid).color = Color::RED();
        gs.m.show_vert_color();
      }
    }
    return false;
  };
}

//=============================== MAIN =========================================

int main(int argc, char **argv) {
  // SETUP GLOBAL STATE AND GUI
  State gs;
  GLcanvas gui = Init_GUI();
  Setup_GUI_Callbacks(gui, gs);
  Setup_Mouse_Callback(gui, gs);

  // Setup font
  ImGui::CreateContext();
  ImGuiIO &io = ImGui::GetIO();
  io.Fonts->Clear(); // Clear any existing fonts
  lato_regular =
      io.Fonts->AddFontFromFileTTF("../font/Lato/Lato-Regular.ttf", 160.0f);
  lato_bold =
      io.Fonts->AddFontFromFileTTF("../font/Lato/Lato-Bold.ttf", 160.0f);
  lato_bold_title =
      io.Fonts->AddFontFromFileTTF("../font/Lato/Lato-Bold.ttf", 180.0f);
  if (lato_regular == NULL || lato_bold == NULL || lato_bold_title == NULL) {
    std::cerr << "Failed to load font" << std::endl;
  }

  // Load mesh
  if (argc > 1) {
    string s = string(argv[1]);
    gs.mesh_path = s;
    Load_mesh(s, gui, gs);
  } else {
    // string s = "../data/pymeshlab_generated/bunny_ok.obj";
    // string s = "../data/cinolib/3holes.obj";
    // string s = "../data/cinolib/bunny.obj";
    // string s = "../pymeshlab/Esperimento_1/data/bob/bob_tri_final.obj";

    // XCODE
    string s = "../../pymeshlab/Esperimento_2/data/15_spot.obj";

    gs.mesh_path = s;
    Load_mesh(s, gui, gs);
  }

  std::thread graph_thread([&]() { init_methods(gs, progress); });
  // do init methods without a separate thread
  // init_methods(gs, progress);

  // Launch the GUI
  gui.launch();

  // Wait for the thread to complete
  graph_thread.join();

  return 0;
}
