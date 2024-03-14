
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

// SSDG with Heat method
#include <cinolib/geodesics.h>

// SSDG with VTP method
#include "SSGD_methods/VTP/vtp_wrapper.h"

// SSGD with graph-based methods
#include "SSGD_methods/Graph-based_methods/extended_solver.h"
using namespace std;
using namespace cinolib;

//------ Global variable ------
// Fonts
ImFont *lato_bold = nullptr;
ImFont *lato_regular = nullptr;
ImFont *lato_bold_title = nullptr;

//:::::::::::::::::::::::::::: GLOBAL VARIABLES (FOR GUI)
//:::::::::::::::::::::::::::::::

struct State {
  // Program state ::::::::::::::::::::::::::::::::::::::::::::::::::::::
  bool MESH_IS_LOADED;
  // Input
  DrawableTrimesh<> m;     // the input mesh
  uint nverts;             // its #vertices
  vector<vector<uint>> VV; // its VV relation

  // GUI state ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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

  // SSGD Method
  enum SSGDMethod { VTP, HEAT, FMM, GEOTANGLE, EDGE } ssgd_method;

  //-------- HEAT method --------
  // sources
  std::vector<uint> sources;
  GeodesicsCache prefactored_matrices;

  //-------- VTP method --------
  vector<double> field_data;
  ScalarField field;
  geodesic_solver solver;
  // sources
  vector<int> voronoi_centers;

  // -------- GeoTangle method --------
  ScalarField field_geo;
  geodesic_solver solver_geo;
  // sources
  vector<int> sources_geo;

  // -------- Edge method --------
  ScalarField field_edge;
  geodesic_solver solver_edge;
  // sources
  vector<int> sources_edge;

  // ----- Timer -----
  double time_heat;
  double vtp_graph_time, vtp_geodesic_time;
  double geotangle_graph_time, geotangle_geodesic_time;
  double edge_graph_time, edge_geodesic_time;

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

    render_mode = State::RENDER_SMOOTH;
    ssgd_method = State::VTP;

    // ----- Timer -----
    time_heat = 0.0;
    vtp_graph_time = vtp_geodesic_time = 0.0;
    geotangle_graph_time = geotangle_geodesic_time = 0.0;
    edge_graph_time = edge_geodesic_time = 0.0;
  }
};

//:::::::::::::::::::::::::::::::::::: I/O ::::::::::::::::::::::::::::::::::::
void Load_mesh(string filename, GLcanvas &gui, State &gs) {
  gs.m = DrawableTrimesh<>(filename.c_str());
  gs.nverts = gs.m.num_verts();
  gs.VV.resize(gs.nverts); // fill in Vertex-Vertex relation
  for (auto i = 0; i < gs.nverts; i++)
    gs.VV[i] = gs.m.vert_ordered_verts_link(i);
  // uncomment the following and adjust parameters if you want a smoother mesh
  // MCF(m,12,1e-5,true);
  gs.m.normalize_bbox(); // rescale mesh to fit [0,1]^3 box
  gs.m.center_bbox();
  gs.m.show_wireframe(gs.SHOW_WIREFRAME);
  gs.m.show_mesh(gs.SHOW_MESH);
  gs.m.updateGL();
  // gs.point_size = gs.m.edge_avg_length()/2; // set initial radius of spheres
  // for critical points

  if (gs.MESH_IS_LOADED) {
    // Clear and reinitialize the vector field for the new mesh
    gs.vec_field = DrawableVectorField();
    gs.show_vecfield = false; // Reset the flag to not show the old vector field
    // Clear the sources for the new mesh
    gs.sources.clear();         // Reset the sources for the new mesh
    gs.voronoi_centers.clear(); // Reset the sources for the new mesh
    gs.sources_geo.clear();     // Reset the sources for the new mesh
    gs.sources_edge.clear();    // Reset the sources for the new mesh
    // Clear cache for Heat method
    gs.prefactored_matrices.heat_flow_cache = NULL; // Reset the heat flow cache

    gs.time_heat = 0.0;
    gs.vtp_graph_time = gs.vtp_geodesic_time = 0.0;
    gs.geotangle_graph_time = gs.geotangle_geodesic_time = 0.0;
    gs.edge_graph_time = gs.edge_geodesic_time = 0.0;
  }

  if (!gs.MESH_IS_LOADED) {
    gui.push(&gs.m);
    gs.MESH_IS_LOADED = true;
  }
}

void Load_mesh(GLcanvas &gui, State &gs) {
  string filename = file_dialog_open();
  if (filename.size() != 0)
    Load_mesh(filename, gui, gs);
}

//::::::::::::::::::::::::::::::::::::: SSGD COMPUTATION
//:::::::::::::::::::::::::::::::::::
void SSGD_Heat(DrawableTrimesh<> &m, GeodesicsCache &prefactored_matrices,
               vector<uint> &sources, double &time_heat) {
  bool cache = false;
  if (prefactored_matrices.heat_flow_cache != NULL) {
    cache = true;
  }
  // Timer
  auto start_heat = chrono::high_resolution_clock::now();
  // Method inside: #include <cinolib/geodesics.h>
  compute_geodesics_amortized(m, prefactored_matrices, sources).copy_to_mesh(m);
  auto stop_heat = chrono::high_resolution_clock::now();
  auto duration_heat =
      chrono::duration_cast<chrono::milliseconds>(stop_heat - start_heat);
  time_heat =
      chrono::duration_cast<chrono::milliseconds>(stop_heat - start_heat)
          .count();
  cout << "Time heat" << time_heat << endl;

  m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);

  if (cache) {
    cout << "Heat computation with cache: " << duration_heat.count()
         << " milliseconds" << endl;
  } else {
    cout << "Heat computation without cache: " << duration_heat.count()
         << " milliseconds" << endl;
  }
}

void SSGD_VTP(DrawableTrimesh<> &m, geodesic_solver &solver,
              vector<double> &field_data, ScalarField &field,
              vector<int> &sources, double &vtp_graph_time,
              double &vtp_geodesic_time) {

  auto start_geodesic_VTP = chrono::high_resolution_clock::now();
  field_data =
      exact_geodesic_distance(m.vector_polys(), m.vector_verts(), sources[0]);
  auto stop_geodesic_VTP = chrono::high_resolution_clock::now();

  // Invert the color mapping
  for (auto &value : field_data) {
    value = 1.0 - value;
  }

  field = ScalarField(field_data);
  field.normalize_in_01();
  field.copy_to_mesh(m);
  m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);

  // auto duration_graph_VTP = chrono::duration_cast<chrono::milliseconds>(
  //     stop_graph_VTP - start_graph_VTP);
  // vtp_graph_time = chrono::duration_cast<chrono::milliseconds>(stop_graph_VTP
  // -
  //                                                              start_graph_VTP)
  //                      .count();
  auto duration_geodesic_VTP = chrono::duration_cast<chrono::milliseconds>(
      stop_geodesic_VTP - start_geodesic_VTP);
  // vtp_geodesic_time = chrono::duration_cast<chrono::milliseconds>(
  //                         stop_geodesic_VTP - start_geodesic_VTP)
  //                         .count();
  // cout << "Graph construction with VTP: " << duration_graph_VTP.count()
  //      << " milliseconds" << endl;
  cout << "Geodesic computation with VTP: " << duration_geodesic_VTP.count()
       << " milliseconds" << endl;
}

// void SSGD_Extended(DrawableTrimesh<> &m, geodesic_solver &solver,
//                    ScalarField &field_geo, vector<int> &sources, const int k,
//                    double &extended_graph_time,
//                    double &extended_geodesic_time) {
//   auto start_graph_extended = chrono::high_resolution_clock::now();

//   solver = extended_solver(m, solver, k);

//   auto stop_graph_extended = chrono::high_resolution_clock::now();

//   vector<double> distances_geo;
//   // type = 0 for geodesic, 1 for isophotic
//   int type = 0;

//   auto start_geodesic_extended = chrono::high_resolution_clock::now();
//   distances_geo = compute_geodesic_distances_geo(solver, sources, type);
//   // update_geodesic_distances_geo(distances_geo, solver, sources, type);
//   auto stop_geodesic_extended = chrono::high_resolution_clock::now();

//   // Invert the color mapping
//   for (auto &value : distances_geo) {
//     value = 1.0 - value;
//   }

//   field_geo = ScalarField(distances_geo);
//   field_geo.normalize_in_01();
//   field_geo.copy_to_mesh(m);
//   m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);

//   auto duration_graph_extended = chrono::duration_cast<chrono::milliseconds>(
//       stop_graph_extended - start_graph_extended);
//   extended_graph_time = chrono::duration_cast<chrono::milliseconds>(
//                             stop_graph_extended - start_graph_extended)
//                             .count();
//   auto duration_geodesic_extended =
//   chrono::duration_cast<chrono::milliseconds>(
//       stop_geodesic_extended - start_geodesic_extended);
//   extended_geodesic_time = chrono::duration_cast<chrono::milliseconds>(
//                                stop_geodesic_extended -
//                                start_geodesic_extended) .count();
//   cout << "Graph construction with extended: "
//        << duration_graph_extended.count() << " milliseconds" << endl;
//   cout << "Geodesic computation with extended: "
//        << duration_geodesic_extended.count() << " milliseconds" << endl;
// }
void SSGD_GeoTangle(DrawableTrimesh<> &m, geodesic_solver &solver,
                    ScalarField &field_geo, vector<int> &sources,
                    double &geotangle_graph_time,
                    double &geotangle_geodesic_time) {
  auto start_graph_GeoTangle = chrono::high_resolution_clock::now();
  solver = make_geodesic_solver(m, true);
  auto stop_graph_GeoTangle = chrono::high_resolution_clock::now();

  vector<double> distances_geo;
  // type = 0 for geodesic, 1 for isophotic
  int type = 0;

  auto start_geodesic_GeoTangle = chrono::high_resolution_clock::now();
  distances_geo = compute_geodesic_distances(solver, sources);
  // update_geodesic_distances_geo(distances_geo, solver, sources, type);
  auto stop_geodesic_GeoTangle = chrono::high_resolution_clock::now();

  // Invert the color mapping
  for (auto &value : distances_geo) {
    value = 1.0 - value;
  }

  field_geo = ScalarField(distances_geo);
  field_geo.normalize_in_01();
  field_geo.copy_to_mesh(m);
  m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);

  auto duration_graph_GeoTangle = chrono::duration_cast<chrono::milliseconds>(
      stop_graph_GeoTangle - start_graph_GeoTangle);
  geotangle_graph_time = chrono::duration_cast<chrono::milliseconds>(
                             stop_graph_GeoTangle - start_graph_GeoTangle)
                             .count();
  auto duration_geodesic_GeoTangle =
      chrono::duration_cast<chrono::milliseconds>(stop_geodesic_GeoTangle -
                                                  start_geodesic_GeoTangle);
  geotangle_geodesic_time =
      chrono::duration_cast<chrono::milliseconds>(stop_geodesic_GeoTangle -
                                                  start_geodesic_GeoTangle)
          .count();
  cout << "Graph construction with GeoTangle: "
       << duration_graph_GeoTangle.count() << " milliseconds" << endl;
  cout << "Geodesic computation with GeoTangle: "
       << duration_geodesic_GeoTangle.count() << " milliseconds" << endl;
}

void SSGD_Edge(DrawableTrimesh<> &m, geodesic_solver &solver,
               ScalarField &field_edge, vector<int> &sources,
               double &edge_graph_time, double &edge_geodesic_time) {
  auto start_graph_edge = chrono::high_resolution_clock::now();
  solver = make_geodesic_solver(m, false);
  auto stop_graph_edge = chrono::high_resolution_clock::now();

  vector<double> distances_edge;
  // type = 0 for geodesic, 1 for isophotic
  int type = 0;

  auto start_geodesic_edge = chrono::high_resolution_clock::now();
  distances_edge = compute_geodesic_distances(solver, sources);
  // update_geodesic_distances_edge(distances_edge, solver, sources, type);
  auto stop_geodesic_edge = chrono::high_resolution_clock::now();

  // Invert the color mapping
  for (auto &value : distances_edge) {
    value = 1.0 - value;
  }

  field_edge = ScalarField(distances_edge);
  field_edge.normalize_in_01();
  field_edge.copy_to_mesh(m);
  m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);

  auto duration_graph_edge = chrono::duration_cast<chrono::milliseconds>(
      stop_graph_edge - start_graph_edge);
  edge_graph_time = chrono::duration_cast<chrono::milliseconds>(
                        stop_graph_edge - start_graph_edge)
                        .count();
  auto duration_geodesic_edge = chrono::duration_cast<chrono::milliseconds>(
      stop_geodesic_edge - start_geodesic_edge);
  edge_geodesic_time = chrono::duration_cast<chrono::milliseconds>(
                           stop_geodesic_edge - start_geodesic_edge)
                           .count();
  cout << "Graph construction with Edge: " << duration_graph_edge.count()
       << " milliseconds" << endl;
  cout << "Geodesic computation with Edge: " << duration_geodesic_edge.count()
       << " milliseconds" << endl;
}

//::::::::::::::::::::::::::::::::::::: GUI
//:::::::::::::::::::::::::::::::::::::::::::::::::
GLcanvas Init_GUI() {
  GLcanvas gui(1500, 700);
  gui.side_bar_width = 0.28;
  gui.show_side_bar = true;
  return gui;
}

void Setup_GUI_Callbacks(GLcanvas &gui, State &gs) {
  gui.callback_app_controls = [&]() {
    // New detached window
    bool show_new_window =
        true; // You can control the visibility with a variable
    if (show_new_window) {
      // Assuming the main window is positioned at (0, 0) and covers the whole
      // screen and considering the sidebar width
      float sidebar_width =
          gui.side_bar_width *
          1500; // Adjust this if side_bar_width is not a ratio
      ImVec2 new_window_pos = ImVec2(
          sidebar_width + 800,
          25); // Top-right corner of the main GUI, right after the sidebar

      ImGui::SetNextWindowPos(new_window_pos, ImGuiCond_FirstUseEver);
      ImVec2 window_size = ImVec2(250, 320); // Example size, change as needed
      ImGui::SetNextWindowSize(window_size, ImGuiCond_FirstUseEver);

      ImGui::Begin("SSGD Methods Timing Results", &show_new_window);

      // Display the label "Sources"
      ImGui::PushFont(lato_bold);
      ImGui::Text("Sources:");
      ImGui::PopFont();
      // Iterate over the gs.sources to display each source vertex ID
      for (uint i = 0; i < gs.sources.size(); ++i) {
        ImGui::Text("Vertex ID: %u", gs.sources[i]);
      }
      ImGui::Text("");

      ImGui::SeparatorText("Exact Polyhedral Methods");
      ImGui::Text("VTP graph time: %.2f ms", gs.vtp_graph_time);
      ImGui::Text("VTP geodesic time: %.2f ms", gs.vtp_geodesic_time);

      ImGui::SeparatorText("PDE-Based Methods");
      ImGui::Text("Heat time: %.2f ms", gs.time_heat);

<<<<<<< Updated upstream
      ImGui::SeparatorText("Graph-Based Methods");
      ImGui::Text("GeoTangle graph time: %.2f ms", gs.geotangle_graph_time);
      ImGui::Text("GeoTangle geodesic time: %.2f ms",
                  gs.geotangle_geodesic_time);
      ImGui::Text("Edge graph time: %.2f ms", gs.edge_graph_time);
      ImGui::Text("Edge geodesic time: %.2f ms", gs.edge_geodesic_time);

      ImGui::End();
=======
        ImGui::End();
>>>>>>> Stashed changes
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
                                 // shaded mesh
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
        // cout << "VTP Method" << endl;
      }
      // PDE-Based Methods
      ImGui::PushFont(lato_bold);
      ImGui::SeparatorText("PDE-Based Methods");
      ImGui::PopFont();
      if (ImGui::RadioButton("Heat  ", gs.ssgd_method == State::HEAT)) {
        gs.ssgd_method = State::HEAT;
        // cout << "HEAT Method" << endl;
      }
      if (ImGui::RadioButton("FMM  ", gs.ssgd_method == State::FMM)) {
        gs.ssgd_method = State::FMM;
        // cout << "FMM Method" << endl;
      }
      // Graph-Based Methods
      ImGui::PushFont(lato_bold);
      ImGui::SeparatorText("Graph-Based Methods");
      ImGui::PopFont();
      if (ImGui::RadioButton("GeoTangle  ",
                             gs.ssgd_method == State::GEOTANGLE)) {
        gs.ssgd_method = State::GEOTANGLE;
        // cout << "GeoTangle Method" << endl;
      }
      if (ImGui::RadioButton("Edge  ", gs.ssgd_method == State::EDGE)) {
        gs.ssgd_method = State::EDGE;
        // cout << "Edge Method" << endl;
      }

      ImGui::TreePop();
    } else {
      ImGui::PopFont();
    }

    ImGui::Text("");
    ImGui::Text("Pick a source vertex: shift + left click");

    // Display the label "Sources"
    ImGui::Text("Sources:");
    // Iterate over the gs.sources to display each source vertex ID
    for (uint i = 0; i < gs.sources.size(); ++i) {
      ImGui::Text("Vertex ID: %u", gs.sources[i]);
    }
    ImGui::Text("");

    // Button for Compute SSGD
    if (ImGui::Button("Compute SSGD")) {
      // Based on the selected SSGD method, perform different actions
      if (gs.sources.empty() && gs.voronoi_centers.empty()) {
        // Open a warning popup if no source is selected
        // cout << "No source selected" << endl;
        ImGui::OpenPopup("Warning");
      } else {

        switch (gs.ssgd_method) {

        case State::VTP: {
          // cout << "Computing SSGD with VTP Method" << endl;
          SSGD_VTP(gs.m, gs.solver, gs.field_data, gs.field, gs.voronoi_centers,
                   gs.vtp_graph_time, gs.vtp_geodesic_time);
          break;
        }

        case State::HEAT: {
          // cout << "Computing SSGD with HEAT Method" << endl;
          SSGD_Heat(gs.m, gs.prefactored_matrices, gs.sources, gs.time_heat);
          break;
        }

        case State::FMM: {
          // cout << "Computing SSGD with FMM Method" << endl;
          //  Add your code for FMM method here
          break;
        }

        case State::GEOTANGLE: {
          // cout << "Computing SSGD with GeoTangle Method" << endl;
          SSGD_GeoTangle(gs.m, gs.solver_geo, gs.field_geo, gs.sources_geo,
                         gs.geotangle_graph_time, gs.geotangle_geodesic_time);
          break;
        }

        case State::EDGE: {
          // cout << "Computing SSGD with Edge Method" << endl;
          SSGD_Edge(gs.m, gs.solver_edge, gs.field_edge, gs.sources_edge,
                    gs.edge_graph_time, gs.edge_geodesic_time);
          break;
        }

        default:
          // cout << "No SSGD method selected" << endl;
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
<<<<<<< Updated upstream
    if (ImGui::SmallButton("Reset")) {
      // Reset HEAT sources
      gs.sources.clear();
      // Reset VTP sources
      gs.voronoi_centers.clear();
      // Reset GeoTangle sources
      gs.sources_geo.clear();
      // Reset Edge sources
      gs.sources_edge.clear();
=======
    if(ImGui::SmallButton("Reset")){
        // Reset HEAT sources
        gs.sources.clear();
        // Reset VTP sources
        gs.voronoi_centers.clear();
        // Reset GeoTangle sources
        gs.sources_geo.clear();
        // Reset Edge sources
        gs.sources_edge.clear();

        gs.time_heat = 0.0;
        gs.vtp_graph_time = gs.vtp_geodesic_time = 0.0;
        gs.geotangle_graph_time = gs.geotangle_geodesic_time = 0.0;
        gs.edge_graph_time = gs.edge_geodesic_time = 0.0;

        // Reset the scalar field
        for(uint vid = 0; vid < gs.m.num_verts(); ++vid) {
          gs.m.vert_data(vid).color = Color::WHITE(); // Replace `original_color` with the actual color
        }
        gs.m.show_vert_color();
>>>>>>> Stashed changes

      gs.time_heat = 0.0;
      gs.vtp_graph_time = 0.0;
      gs.vtp_geodesic_time = 0.0;
      gs.geotangle_graph_time = 0.0;
      gs.geotangle_geodesic_time = 0.0;
      gs.edge_graph_time = 0.0;
      gs.edge_geodesic_time = 0.0;

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
<<<<<<< Updated upstream
  gui.callback_mouse_left_click = [&](int modifiers) -> bool {
    if (modifiers & GLFW_MOD_SHIFT) {
      vec3d p;
      vec2d click = gui.cursor_pos();
      if (gui.unproject(click, p)) {
        //----- HEAT SOURCES -----
        uint vid = gs.m.pick_vert(p);
        gs.sources.push_back(vid);
        cout << "Source vertex: " << vid << endl;
=======
    gui.callback_mouse_left_click = [&](int modifiers) -> bool {
        if(modifiers & GLFW_MOD_SHIFT) {
            vec3d p;
            vec2d click = gui.cursor_pos();
            if(gui.unproject(click, p)) {
                // ----- HEAT SOURCES -----
                uint vid = gs.m.pick_vert(p);
                gs.sources.push_back(vid);
                cout << "Source vertex: " << vid << endl;
>>>>>>> Stashed changes

        // Color the selected vertex in red
        gs.m.vert_data(vid).color = Color::RED();
        gs.m.show_vert_color();

<<<<<<< Updated upstream
        //------ VTP SOURCES ------
        int selected_vid = gs.m.pick_vert(p);
        gs.voronoi_centers.push_back(selected_vid);
        // std::cout << "Selected vid VTP = " << selected_vid << std::endl;

        // GeoTangle sources
        gs.sources_geo.push_back(selected_vid);
        // std::cout << "Selected vid GEO = " << selected_vid << std::endl;

        // Edge sources
        gs.sources_edge.push_back(selected_vid);
        // std::cout << "Selected vid EDGE = " << selected_vid << std::endl;

        // You might need to replace "profiler" with your own profiling method
        // or remove it profiler.push("compute_geodesics");
        // compute_geodesics_amortized(gs.m, gs.prefactored_matrices,
        // gs.sources).copy_to_mesh(gs.m);
        // profiler.pop();
        // gs.m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);
      }
    }
    return false;
  };
=======
                // ------ VTP SOURCES ------
                int selected_vid = gs.m.pick_vert(p);
                gs.voronoi_centers.push_back(selected_vid);
                //std::cout << "Selected vid VTP = " << selected_vid << std::endl;

                // ------ GEOTANGLE SOURCES ------
                gs.sources_geo.push_back(selected_vid);
                //std::cout << "Selected vid GEO = " << selected_vid << std::endl;

                // ------ EDGE SOURCES ------
                gs.sources_edge.push_back(selected_vid);
                //std::cout << "Selected vid EDGE = " << selected_vid << std::endl;
            }
        }
        return false;
    };
>>>>>>> Stashed changes
}

//=============================== MAIN =========================================

int main(int argc, char **argv) {

  // SETUP GLOBAL STATE AND GUI:::::::::::::::::::::::
  State gs;
  GLcanvas gui = Init_GUI();
  Setup_GUI_Callbacks(gui, gs);
  Setup_Mouse_Callback(gui, gs);

  // Setup font
  ImGui::CreateContext();
  ImGuiIO &io = ImGui::GetIO();
  io.Fonts->Clear(); // Clear any existing fonts
<<<<<<< Updated upstream
  lato_regular =
      io.Fonts->AddFontFromFileTTF("../font/Lato/Lato-Regular.ttf", 160.0f);
  lato_bold =
      io.Fonts->AddFontFromFileTTF("../font/Lato/Lato-Bold.ttf", 160.0f);
  lato_bold_title =
      io.Fonts->AddFontFromFileTTF("../font/Lato/Lato-Bold.ttf", 180.0f);
  if (lato_regular == NULL || lato_bold == NULL) {
    std::cerr << "Failed to load font" << std::endl;
  }

  // Load mesh
  if (argc > 1) {
=======
  lato_regular = io.Fonts->AddFontFromFileTTF("../font/Lato/Lato-Regular.ttf", 160.0f);
  lato_bold = io.Fonts->AddFontFromFileTTF("../font/Lato/Lato-Bold.ttf", 160.0f);
  lato_bold_title = io.Fonts->AddFontFromFileTTF("../font/Lato/Lato-Bold.ttf", 180.0f);
  if(lato_regular == NULL || lato_bold == NULL || lato_bold_title == NULL) {
      std::cerr << "Failed to load font" << std::endl;
  }

  //Load mesh
  if (argc>1) {
>>>>>>> Stashed changes
    string s = "../data/" + string(argv[1]);
    Load_mesh(s, gui, gs);
  } else {
    string s = "../data/cinolib/bunny.obj";
    Load_mesh(s, gui, gs);
  }

  // // // GENERATE FIELD:::::::::::::::::::::::::::::::
  // Generate_field(gui,gs);

  // // COMPUTE DISCRETE SCALE SPACE:::::::::::::::::
  // Build_disc_ss(gui,gs);

  // render the mesh
  return gui.launch();
}
