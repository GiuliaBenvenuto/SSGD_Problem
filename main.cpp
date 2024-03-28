
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
using namespace std;
using namespace cinolib;

// Compute SSGD
#include "solving_ssgd.h"

//------ Global variable ------
// Fonts
ImFont *lato_bold = nullptr;
ImFont *lato_regular = nullptr;
ImFont *lato_bold_title = nullptr;
atomic<float> progress(0.0f);
bool isSettingK = false;
atomic<bool> isSolverThreadRunning(false);
std::future<void> solver_future;



//:::::::::::::::::::::::::::: GLOBAL VARIABLES (FOR GUI):::::::::::::::::::::::::::::::
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

  //-------- GUI state --------
  // View
  bool SHOW_MESH, SHOW_WIREFRAME;

  // Wireframe
  float wireframe_width; // Wireframe width
  float wireframe_alpha; // Wireframe transparency

  // Shading
  enum MeshRenderMode {RENDER_POINTS, RENDER_FLAT, RENDER_SMOOTH} render_mode;

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
  enum SSGDMethod {VTP, TRETTNER, HEAT, FMM, GEOTANGLE, EDGE, EXTENDED} ssgd_method;
  
  ScalarField field;

  // Cache for Heat method
  GeodesicsCache prefactored_matrices;

  // K for extended
  int k;
  int prev_k;

  // Time for heat
  float heat_time;
  float heat_time_prev;

  // Sources for SSGD
  std::vector<uint> sources_heat;
  vector<int> sources;

  // Trettner
  string mesh_path;
  HalfEdge mesh;

  VTPSolver         vtp_solver;
  TrettnerSolver    trettner_solver; 
  HeatSolver        heat_solver;
  GeotangleSolver   geotangle_solver;
  EdgeSolver        edge_solver;
  ExtendedSolver    extended_solver;
  

  // Timer
  std::chrono::steady_clock::time_point tic;
  std::chrono::steady_clock::time_point toc;
  double vtp_load,        vtp_preprocess,         vtp_query;
  double trettner_load,   trettner_preprocess,    trettner_query;
  double heat_load,       heat_preprocess,        heat_query;
  double geotangle_load,  geotangle_preprocess,   geotangle_query;
  double edge_load,       edge_preprocess,        edge_query;
  double extended_load,   extended_preprocess,    extended_query;


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
    vtp_load,       vtp_preprocess,       vtp_query = 0.0;
    trettner_load,  trettner_preprocess,  trettner_query = 0.0;
    heat_load,      heat_preprocess,      heat_query = 0.0;
    geotangle_load, geotangle_preprocess, geotangle_query = 0.0;
    edge_load,      edge_preprocess,       edge_query = 0.0;
    extended_load,  extended_preprocess,  extended_query = 0.0;

    res = vector<double>();

    mesh_path = "";
    mesh = HalfEdge();

    k = 3;
    prev_k = 3;
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

void fillTimeTable(State &gs, const string& method_name, double load_time, double preprocess_time, double query_time) {
  if (method_name == "VTP") {
    gs.vtp_load = load_time;
    gs.vtp_preprocess = preprocess_time;
    gs.vtp_query = query_time;
  } else if (method_name == "Trettner") {
    gs.trettner_load = load_time;
    gs.trettner_preprocess = preprocess_time;
    gs.trettner_query = query_time;
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
  }
}


void init(GeodesicMethod &m, State &gs, string name) {
  cout << "---------- Initializing method: " << name << " ----------" << endl;

  // load
  gs.tic = std::chrono::steady_clock::now();
  m.load(gs.coords, gs.tris);
  gs.toc = std::chrono::steady_clock::now();
  double load = chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic).count();
  cout << "Load time: " << load << " milliseconds" << endl;

  // preprocess
  gs.tic = std::chrono::steady_clock::now();
  m.preprocess();
  gs.toc = std::chrono::steady_clock::now();
  double preprocess = chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic).count();
  cout << "Preprocess time: " << preprocess << " milliseconds" << endl;

  // fill time table
  fillTimeTable(gs, name, load, preprocess, 0.0);
}


void init_methods(State &gs, atomic<float> &progress) {
  init(gs.vtp_solver, gs, "VTP");
  progress.store(1.0f / 6.0f);  

  init(gs.trettner_solver, gs, "Trettner");
  progress.store(2.0f / 6.0f); 

  init(gs.heat_solver, gs, "Heat");
  progress.store(3.0f / 6.0f);  

  init(gs.geotangle_solver, gs, "Geotangle");
  progress.store(4.0f / 6.0f);  

  init(gs.edge_solver, gs, "Edge");
  progress.store(5.0f / 6.0f);  

  init(gs.extended_solver, gs, "Extended");
  progress.store(1.0f);         
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

    gs.vtp_load,        gs.vtp_preprocess,        gs.vtp_query = 0.0;
    gs.trettner_load,   gs.trettner_preprocess,   gs.trettner_query = 0.0;
    gs.heat_load,       gs.heat_preprocess,       gs.heat_query = 0.0;
    gs.geotangle_load,  gs.geotangle_preprocess,  gs.geotangle_query = 0.0;
    gs.edge_load,       gs.edge_preprocess,       gs.edge_query = 0.0;
    gs.extended_load,   gs.extended_preprocess,   gs.extended_query = 0.0;
    gs.res = vector<double>();

    gs.mesh_path = filename;
    gs.trettner_solver = TrettnerSolver(gs.mesh_path);

    init_methods(gs, progress);
  }

  if (!gs.MESH_IS_LOADED) {
    gs.mesh_path = filename;
    gs.trettner_solver = TrettnerSolver(gs.mesh_path);

    gui.push(&gs.m);
    gs.MESH_IS_LOADED = true;
  }
}

void Load_mesh(GLcanvas &gui, State &gs) {
  string filename = file_dialog_open();
  if (filename.size() != 0)
    gs.mesh_path = filename;
  Load_mesh(filename, gui, gs);
}

//:::::::::::::::::::::::::::::::::::::: GUI:::::::::::::::::::::::::::::::::::::::::::::::::
GLcanvas Init_GUI() {
  GLcanvas gui(1500, 700);
  gui.side_bar_width = 0.28;
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
          ImGui::Text("%.2f", gs.vtp_load);
          ImGui::TableSetColumnIndex(2);
          ImGui::Text("%.2f", gs.vtp_preprocess);
          ImGui::TableSetColumnIndex(3);
          ImGui::Text("%.2f", gs.vtp_query);

          // Trettner Method
          ImGui::TableNextRow();
          ImGui::TableSetColumnIndex(0);
          ImGui::Text("Trettner");
          ImGui::TableSetColumnIndex(1);
          ImGui::Text("%.2f", gs.trettner_load);
          ImGui::TableSetColumnIndex(2);
          ImGui::Text("%.2f", gs.trettner_preprocess);
          ImGui::TableSetColumnIndex(3);
          ImGui::Text("%.2f", gs.trettner_query);

          // Heat Method
          ImGui::TableNextRow();
          ImGui::TableSetColumnIndex(0);
          ImGui::Text("Heat");
          ImGui::TableSetColumnIndex(1);
          ImGui::Text("%.2f", gs.heat_load);
          ImGui::TableSetColumnIndex(2);
          ImGui::Text("%.2f", gs.heat_preprocess);
          ImGui::TableSetColumnIndex(3);
          ImGui::Text("%.2f", gs.heat_query);

          // GeoTangle Method
          ImGui::TableNextRow();
          ImGui::TableSetColumnIndex(0);
          ImGui::Text("GeoTangle");
          ImGui::TableSetColumnIndex(1);
          ImGui::Text("%.2f", gs.geotangle_load);
          ImGui::TableSetColumnIndex(2);
          ImGui::Text("%.2f", gs.geotangle_preprocess);
          ImGui::TableSetColumnIndex(3);
          ImGui::Text("%.2f", gs.geotangle_query);

          // Edge Method
          ImGui::TableNextRow();
          ImGui::TableSetColumnIndex(0);
          ImGui::Text("Edge");
          ImGui::TableSetColumnIndex(1);
          ImGui::Text("%.2f", gs.edge_load);
          ImGui::TableSetColumnIndex(2);
          ImGui::Text("%.2f", gs.edge_preprocess);
          ImGui::TableSetColumnIndex(3);
          ImGui::Text("%.2f", gs.edge_query);

          // Extended Method
          ImGui::TableNextRow();
          ImGui::TableSetColumnIndex(0);
          ImGui::Text("Extended");
          ImGui::TableSetColumnIndex(1);
          ImGui::Text("%.2f", gs.extended_load);
          ImGui::TableSetColumnIndex(2);
          ImGui::Text("%.2f", gs.extended_preprocess);
          ImGui::TableSetColumnIndex(3);
          ImGui::Text("%.2f", gs.extended_query);

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
      ImGui::ProgressBar(progress.load(), ImVec2(-1.0f, 0.0f), progress > 1.0f ? "Done" : "Loading...");
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
      int numFaces = gs.m.num_polys(); // or num_faces() depending on your mesh type

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
      // if (ImGui::RadioButton("Heat  ", gs.ssgd_method == State::HEAT)) {
      //   gs.ssgd_method = State::HEAT;
      // }
      // ImGui::SameLine(0, 60);
      // ImGui::SetNextItemWidth(gs.width);
      // ImGui::InputFloat("param", &gs.heat_time, 0.1f, 1.0f, "%.3f");
      // Define two columns
      ImGui::Columns(2, nullptr, false);
      // First column for the radio button
      if (ImGui::RadioButton("Heat  ", gs.ssgd_method == State::HEAT)) {
          gs.ssgd_method = State::HEAT;
      }
      ImGui::NextColumn();
      // Second column for the InputFloat
      ImGui::SetNextItemWidth(gs.width); 
      ImGui::InputFloat("t", &gs.heat_time, 0.1f, 1.0f, "%.3f");
      ImGui::Columns(1);

      if (ImGui::RadioButton("FMM  ", gs.ssgd_method == State::FMM)) {
        gs.ssgd_method = State::FMM;
      }

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
      // Second column for the InputFloat
      ImGui::SetNextItemWidth(gs.width); 
      ImGui::InputInt("k", &gs.k, 1, 10); // Add an input integer for parameter k
      ImGui::Columns(1);
      gs.k = std::max(gs.k, 1);

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
          gs.vtp_solver.query(gs.sources[0], gs.res, gs.field);
          gs.toc = std::chrono::steady_clock::now();
          gs.vtp_query = chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic).count();
          fillTimeTable(gs, "VTP", gs.vtp_load, gs.vtp_preprocess, gs.vtp_query);

          gs.field.copy_to_mesh(gs.m);
          gs.m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);
          break;
        }

        case State::TRETTNER: {
          gs.tic = std::chrono::steady_clock::now();
          gs.trettner_solver.query(gs.sources[0], gs.res, gs.field);
          gs.toc = std::chrono::steady_clock::now();
          gs.trettner_query = chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic).count();
          fillTimeTable(gs, "Trettner", gs.trettner_load, gs.trettner_preprocess, gs.trettner_query);

          gs.field.copy_to_mesh(gs.m);
          gs.m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);
          break;
        }

        case State::HEAT: {
          if (gs.heat_time != gs.heat_time_prev) {
            cout << "Heat time has changed to: " << gs.heat_time << endl;
            gs.heat_solver.set_t(gs.heat_time);
            gs.heat_time_prev = gs.heat_time;
          }

          gs.tic = std::chrono::steady_clock::now();
          gs.heat_solver.query(gs.sources[0], gs.res, gs.field);
          gs.toc = std::chrono::steady_clock::now();
          gs.heat_query = chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic).count();
          fillTimeTable(gs, "Heat", gs.heat_load, gs.heat_preprocess, gs.heat_query);

          gs.field.copy_to_mesh(gs.m);
          gs.m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);
          break;
        }

        case State::FMM: {
          // cout << "Computing SSGD with FMM Method" << endl;
          break;
        }

        case State::GEOTANGLE: {
          gs.tic = std::chrono::steady_clock::now();
          gs.geotangle_solver.query(gs.sources[0], gs.res, gs.field);
          gs.toc = std::chrono::steady_clock::now();
          gs.geotangle_query = chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic).count();
          fillTimeTable(gs, "Geotangle", gs.geotangle_load, gs.geotangle_preprocess, gs.geotangle_query);

          gs.field.copy_to_mesh(gs.m);
          gs.m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);
          break;
        }

        case State::EDGE: {
          gs.tic = std::chrono::steady_clock::now();
          gs.edge_solver.query(gs.sources[0], gs.res, gs.field);
          gs.toc = std::chrono::steady_clock::now();
          gs.edge_query = chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic).count();
          fillTimeTable(gs, "Edge", gs.edge_load, gs.edge_preprocess, gs.edge_query);

          gs.field.copy_to_mesh(gs.m);
          gs.m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);
          break;
        }

        case State::EXTENDED: {
          if (gs.k != gs.prev_k) {
            cout << "k has changed to: " << gs.k << endl;
            isSettingK = true;
            
            // Start a new thread to update the solver
            std::thread solver_thread([&gs, isSolverThreadRunningPtr = &isSolverThreadRunning]() {
                gs.tic = std::chrono::steady_clock::now();
                gs.extended_solver.set_k(gs.k); // This is the blocking call
                gs.toc = std::chrono::steady_clock::now();

                // Calculate the time taken for preprocessing
                gs.extended_preprocess = chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic).count();
                fillTimeTable(gs, "Extended", gs.extended_load, gs.extended_preprocess, gs.extended_query);

                gs.prev_k = gs.k;
                isSettingK = false;
                *isSolverThreadRunningPtr = false;
              
            });
            solver_thread.join();
          }

          gs.tic = std::chrono::steady_clock::now();
          gs.extended_solver.query(gs.sources[0], gs.res, gs.field);
          gs.toc = std::chrono::steady_clock::now();
          gs.extended_query = chrono::duration_cast<chrono::milliseconds>(gs.toc - gs.tic).count();
          fillTimeTable(gs, "Extended", gs.extended_load, gs.extended_preprocess, gs.extended_query);

          gs.field.copy_to_mesh(gs.m);
          gs.m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);

          break;
        }
        
        default:
          cout << "No SSGD method selected" << endl;
          break;
        }
      }
    }



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

      gs.vtp_query = 0.0;
      gs.trettner_query = 0.0;
      gs.heat_query = 0.0;
      gs.geotangle_query = 0.0;
      gs.edge_query = 0.0;
      gs.extended_query = 0.0;

      gs.res = vector<double>();

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
    string s = "../data/" + string(argv[1]);
    Load_mesh(s, gui, gs);
  } else {
    string s = "../data/cinolib/bunny.obj";
    // string s = "../data/Trettner/69930.obj";
    gs.mesh_path = s;
    Load_mesh(s, gui, gs);
  }

  std::thread graph_thread([&]() { init_methods(gs, progress); });

  // Launch the GUI
  gui.launch();

  // Wait for the thread to complete
  graph_thread.join();

  return 0;
}
