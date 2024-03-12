// #include "utilities.h"
#include <cinolib/drawable_vector_field.h>
#include <cinolib/drawable_segment_soup.h>
#include <cinolib/drawable_sphere.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/file_dialog_open.h>
#include <cinolib/gl/file_dialog_save.h>
#include <cinolib/gradient.h>
#include <cinolib/scalar_field.h>
#include <cinolib/meshes/drawable_tetmesh.h>
#include <cinolib/vector_serialization.h>
#include <cinolib/io/write_OBJ.h>
#include <cinolib/mean_curv_flow.h>
#include <Eigen/SparseCholesky>
#include <fstream>
#include <imgui.h>
#include <chrono>

// SSDG with Heat method
#include <cinolib/geodesics.h>

// SSDG with VTP method
#include "SSGD_methods/VTP/diff_geo.h"
#include "SSGD_methods/VTP/diff_geo.cpp"

// SSGD with GeoTangle 
#include "SSGD_methods/GeoTangle/GeoTangle.cpp"

// SSGD with Edge method
#include "SSGD_methods/Edge/edge.cpp"

using namespace std;
using namespace cinolib;

//------ Global variable ------
// Fonts
ImFont* lato_bold = nullptr;
ImFont* lato_regular = nullptr;
ImFont* lato_bold_title = nullptr;


//:::::::::::::::::::::::::::: GLOBAL VARIABLES (FOR GUI) ::::::::::::::::::::::::::::::

struct State {
  // Program state ::::::::::::::::::::::::::::::::::::::::::::::::::::::
  bool MESH_IS_LOADED;
  // Input
  DrawableTrimesh<> m;            // the input mesh
  uint nverts;                    // its #vertices
  vector<vector<uint>> VV;        // its VV relation

  // GUI state ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  // View
  bool SHOW_MESH, SHOW_WIREFRAME;
  
  // Wireframe
  float wireframe_width;  // Wireframe width
  float wireframe_alpha;  // Wireframe transparency
  
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
  vector<string> metric_names = {"Geodesic", "Isophotic"};

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



  State() {
    MESH_IS_LOADED = false;
    // view
    SHOW_WIREFRAME = false;
    SHOW_MESH = true;

    wireframe_width = 1.0;
    wireframe_alpha = 1.0;

    // vector field
    show_vecfield = false;
    vecfield_size = 0.9f;
    vec_color = Color::RED();

    vert_color = Color::WHITE();
    poly_color = Color::WHITE();

    render_mode = State::RENDER_SMOOTH;
    ssgd_method = State::VTP;
  }
};



//:::::::::::::::::::::::::::::::::::: I/O ::::::::::::::::::::::::::::::::::::
void Load_mesh(string filename, GLcanvas & gui, State &gs)
{
  gs.m = DrawableTrimesh<>(filename.c_str());
  gs.nverts = gs.m.num_verts();
  gs.VV.resize(gs.nverts); // fill in Vertex-Vertex relation
  for (auto i=0;i<gs.nverts;i++) gs.VV[i]=gs.m.vert_ordered_verts_link(i);
  // uncomment the following and adjust parameters if you want a smoother mesh
  // MCF(m,12,1e-5,true);
  gs.m.normalize_bbox(); // rescale mesh to fit [0,1]^3 box
  gs.m.center_bbox();
  gs.m.show_wireframe(gs.SHOW_WIREFRAME);
  gs.m.show_mesh(gs.SHOW_MESH);    
  gs.m.updateGL();  
  // gs.point_size = gs.m.edge_avg_length()/2; // set initial radius of spheres for critical points

  if (gs.MESH_IS_LOADED) {
    // Clear and reinitialize the vector field for the new mesh
    gs.vec_field = DrawableVectorField();
    gs.show_vecfield = false; // Reset the flag to not show the old vector field
  }

  if (!gs.MESH_IS_LOADED) {
    gui.push(&gs.m);
    gs.MESH_IS_LOADED = true;
  }
}

void Load_mesh(GLcanvas & gui, State &gs)
{
  string filename = file_dialog_open();
  if (filename.size()!=0) Load_mesh(filename,gui,gs);
}


//::::::::::::::::::::::::::::::::::::: SSGD COMPUTATION ::::::::::::::::::::::::::::::::::
void SSGD_Heat(DrawableTrimesh<> &m, GeodesicsCache &prefactored_matrices, vector<uint> &sources) {
  bool cache = false;
  if (prefactored_matrices.heat_flow_cache != NULL) {
    cache = true;
  }
  // Time
  auto start = chrono::high_resolution_clock::now();
  // Method inside: #include <cinolib/geodesics.h>
  compute_geodesics_amortized(m, prefactored_matrices, sources).copy_to_mesh(m);
  auto stop = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);

  m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);
  
  if (cache) {
    cout << "Heat computation with cache: " << duration.count() << " milliseconds" << endl;
  } else {
    cout << "Heat computation without cache: " << duration.count() << " milliseconds" << endl;
  }
}


void SSGD_VTP(DrawableTrimesh<> &m, geodesic_solver &solver, vector<double> &field_data, ScalarField &field, vector<int> &sources) {
  vector<patch> quadrics = patch_fitting(m, 5);
  // Method inside: #include "SSGD_methods/VTP/diff_geo.cpp"
  solver = compute_geodesic_solver(m, quadrics);
  // Geodesic = 0, Isophotic = 1
  const int type_of_metric = 0;
  field_data = compute_geodesic_distances(solver, sources, type_of_metric);

  // Invert the color mapping
  for (auto& value : field_data) {
    value = 1.0 - value;
  }

  field = ScalarField(field_data);
  field.normalize_in_01();
  field.copy_to_mesh(m);
  m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);
}


void SSGD_GeoTangle(DrawableTrimesh<> &m, geodesic_solver &solver, ScalarField &field_geo, vector<int> &sources) {
  solver = make_geodesic_solver(m);
  vector<double> distances_geo;
  // type = 0 for geodesic, 1 for isophotic
  int type = 0;
  distances_geo = compute_geodesic_distances_geo(solver, sources, type);
  update_geodesic_distances_geo(distances_geo, solver, sources, type);

  // Invert the color mapping
  for (auto& value : distances_geo) {
    value = 1.0 - value;
  }

  field_geo = ScalarField(distances_geo);
  field_geo.normalize_in_01();
  field_geo.copy_to_mesh(m);
  m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);
}


void SSGD_Edge(DrawableTrimesh<> &m, geodesic_solver &solver, ScalarField &field_edge, vector<int> &sources) {
  solver = make_geodesic_solver_edge(m);
  vector<double> distances_edge;
  // type = 0 for geodesic, 1 for isophotic
  int type = 0;
  distances_edge = compute_geodesic_distances_edge(solver, sources, type);
  cout << "Distances edge size: " << distances_edge.size() << endl;
  update_geodesic_distances_edge(distances_edge, solver, sources, type);

  // Invert the color mapping
  for (auto& value : distances_edge) {
    value = 1.0 - value;
  }

  field_edge = ScalarField(distances_edge);
  field_edge.normalize_in_01();
  field_edge.copy_to_mesh(m);
  m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);
}



//::::::::::::::::::::::::::::::::::::: GUI ::::::::::::::::::::::::::::::::::::::::::::::::
GLcanvas Init_GUI()
{
  GLcanvas gui(1500,700);
  gui.side_bar_width = 0.28;
  gui.show_side_bar = true;
  return gui;
}

void Setup_GUI_Callbacks(GLcanvas & gui, State &gs)
{
  gui.callback_app_controls = [&]() {
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
          Load_mesh(gui,gs);
        }
      }
      ImGui::SameLine();
      if(ImGui::SmallButton("Save"))
      {
          std::string filename = file_dialog_save();
          if(!filename.empty())
          {
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
      ImGui::SetNextWindowPos(center,ImGuiCond_Appearing,ImVec2(0.5f,0.5f));
      if (ImGui::BeginPopupModal("Load mesh?", NULL, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoCollapse))
      {
        static bool dont_ask_me_next_time = false;
        if (dont_ask_me_next_time) {Load_mesh(gui,gs); ImGui::CloseCurrentPopup();}
        ImGui::Text("All data structures will be reset - Load anyway?\n\n");
        ImGui::Separator();           
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(0,0));
        ImGui::Checkbox("Don't ask me next time", &dont_ask_me_next_time);
        ImGui::PopStyleVar();
        if (ImGui::Button("OK", ImVec2(120,0))) {Load_mesh(gui,gs); ImGui::CloseCurrentPopup();}
        ImGui::SameLine();
        if (ImGui::Button("Cancel", ImVec2(120,0))) {ImGui::CloseCurrentPopup();}
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
      if (ImGui::Checkbox("Show", &gs.SHOW_WIREFRAME)) gs.m.show_wireframe(gs.SHOW_WIREFRAME);
      if (ImGui::SliderFloat("Width", &gs.wireframe_width, 1.f, 10.f)) gs.m.show_wireframe_width(gs.wireframe_width);
      if (ImGui::SliderFloat("Transparency", &gs.wireframe_alpha, 0.f, 1.f)) gs.m.show_wireframe_transparency(gs.wireframe_alpha);
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
      if(ImGui::RadioButton("Point ", gs.render_mode == State::RENDER_POINTS)) {
        gs.render_mode = State::RENDER_POINTS;
        gs.m.show_mesh_points();
      }
      // Flat shading mode
      if(ImGui::RadioButton("Flat  ", gs.render_mode == State::RENDER_FLAT)) {
        gs.render_mode = State::RENDER_FLAT;
        gs.m.show_mesh_flat(); 
      }
      // Smooth shading mode
      if(ImGui::RadioButton("Smooth", gs.render_mode == State::RENDER_SMOOTH)) {
        gs.render_mode = State::RENDER_SMOOTH;
        gs.m.show_mesh_smooth(); // Assuming you have a method to display smooth shaded mesh
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
            if (ImGui::RadioButton("Vert", gs.m.drawlist.draw_mode & DRAW_TRI_VERTCOLOR)) {
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
            if (ImGui::RadioButton("Poly", gs.m.drawlist.draw_mode & DRAW_TRI_FACECOLOR)) {
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
                if (gs.vec_field.size() == 0) { // Initialize the vector field if not already
                    gs.vec_field = DrawableVectorField(gs.m);
                    ScalarField f(gs.m.serialize_uvw(3)); // Adjust U_param if needed
                    gs.vec_field = gradient_matrix(gs.m) * f;
                    gs.vec_field.normalize();
                    // print the values inside vec_field
                    // Print the number of elements in the vector field
                    gs.vec_field.set_arrow_size(float(gs.m.edge_avg_length())*gs.vecfield_size);
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
            gs.vec_field.set_arrow_size(float(gs.m.edge_avg_length()) * gs.vecfield_size);
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
      if(ImGui::RadioButton("VTP ", gs.ssgd_method == State::VTP)) {
        gs.ssgd_method = State::VTP;
        //cout << "VTP Method" << endl;
      }
      // PDE-Based Methods
      ImGui::PushFont(lato_bold);
      ImGui::SeparatorText("PDE-Based Methods");
      ImGui::PopFont();
      if(ImGui::RadioButton("Heat  ", gs.ssgd_method == State::HEAT)) {
        gs.ssgd_method = State::HEAT;
        //cout << "HEAT Method" << endl;
      }
      if(ImGui::RadioButton("FMM  ", gs.ssgd_method == State::FMM)) {
        gs.ssgd_method = State::FMM;
        //cout << "FMM Method" << endl;
      }
      // Graph-Based Methods
      ImGui::PushFont(lato_bold);
      ImGui::SeparatorText("Graph-Based Methods");
      ImGui::PopFont();
      if(ImGui::RadioButton("GeoTangle  ", gs.ssgd_method == State::GEOTANGLE)) {
        gs.ssgd_method = State::GEOTANGLE;
        //cout << "GeoTangle Method" << endl;
      }
      if(ImGui::RadioButton("Edge  ", gs.ssgd_method == State::EDGE)) {
        gs.ssgd_method = State::EDGE;
        //cout << "Edge Method" << endl;
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
            //cout << "No source selected" << endl;
            ImGui::OpenPopup("Warning");
        } else {

          switch (gs.ssgd_method) {
      
            case State::VTP: {
              //cout << "Computing SSGD with VTP Method" << endl;
              SSGD_VTP(gs.m, gs.solver, gs.field_data, gs.field, gs.voronoi_centers);
              break;
            }

            case State::HEAT: {
              //cout << "Computing SSGD with HEAT Method" << endl;
              SSGD_Heat(gs.m, gs.prefactored_matrices, gs.sources);
              break;
            }

            case State::FMM: {
              //cout << "Computing SSGD with FMM Method" << endl;
              // Add your code for FMM method here
              break;
            }

            case State::GEOTANGLE: {
              //cout << "Computing SSGD with GeoTangle Method" << endl;
              SSGD_GeoTangle(gs.m, gs.solver_geo, gs.field_geo, gs.sources_geo);
              break;
            }

            case State::EDGE: {
              //cout << "Computing SSGD with Edge Method" << endl;
              SSGD_Edge(gs.m, gs.solver_edge, gs.field_edge, gs.sources_edge);
              break;
            }

            default:
              //cout << "No SSGD method selected" << endl;
              break;
          }
      }
    }

    // Popup modal for warning if no source is selected
    if (ImGui::BeginPopupModal("Warning", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::Text("Select a source vertex before computing SSGD");
        ImGui::Text("shift + left click");
        if (ImGui::Button("OK")) {
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }

    // Reset button
    if(ImGui::SmallButton("Reset")){
        // Reset HEAT sources
        gs.sources.clear();
        // Reset VTP sources
        gs.voronoi_centers.clear();
        // Reset GeoTangle sources
        gs.sources_geo.clear();
        // Reset Edge sources
        gs.sources_edge.clear();

        // Reset the scalar field
        for(uint vid = 0; vid < gs.m.num_verts(); ++vid) {
          gs.m.vert_data(vid).color = Color::WHITE(); // Replace `original_color` with the actual color
        }
        gs.m.show_vert_color();
    }

  };
}

// Compute the geodesic distances from a single source vertex
void Setup_Mouse_Callback(GLcanvas &gui, State &gs) {
    gui.callback_mouse_left_click = [&](int modifiers) -> bool {
        if(modifiers & GLFW_MOD_SHIFT) {
            vec3d p;
            vec2d click = gui.cursor_pos();
            if(gui.unproject(click, p)) {
                //----- HEAT SOURCES -----
                uint vid = gs.m.pick_vert(p);
                gs.sources.push_back(vid);
                cout << "Source vertex: " << vid << endl;

                // Color the selected vertex in red
                gs.m.vert_data(vid).color = Color::RED();
                gs.m.show_vert_color();

                //------ VTP SOURCES ------
                int selected_vid = gs.m.pick_vert(p);
                gs.voronoi_centers.push_back(selected_vid);
                //std::cout << "Selected vid VTP = " << selected_vid << std::endl;

                // GeoTangle sources
                gs.sources_geo.push_back(selected_vid);
                //std::cout << "Selected vid GEO = " << selected_vid << std::endl;

                // Edge sources
                gs.sources_edge.push_back(selected_vid);
                //std::cout << "Selected vid EDGE = " << selected_vid << std::endl;


                // You might need to replace "profiler" with your own profiling method or remove it
                // profiler.push("compute_geodesics");
                //compute_geodesics_amortized(gs.m, gs.prefactored_matrices, gs.sources).copy_to_mesh(gs.m);
                // profiler.pop();
                //gs.m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);
            }
        }
        return false;
    };
}


//=============================== MAIN =========================================

int main(int argc, char **argv) {

  //SETUP GLOBAL STATE AND GUI:::::::::::::::::::::::
  State gs;
  GLcanvas gui = Init_GUI();
  Setup_GUI_Callbacks(gui, gs);
  Setup_Mouse_Callback(gui, gs);

  // Setup font
  ImGui::CreateContext();
  ImGuiIO& io = ImGui::GetIO();
  io.Fonts->Clear(); // Clear any existing fonts
  lato_regular = io.Fonts->AddFontFromFileTTF("../font/Lato/Lato-Regular.ttf", 160.0f);
  lato_bold = io.Fonts->AddFontFromFileTTF("../font/Lato/Lato-Bold.ttf", 160.0f);
  lato_bold_title = io.Fonts->AddFontFromFileTTF("../font/Lato/Lato-Bold.ttf", 180.0f);
  if(lato_regular == NULL || lato_bold == NULL) {
      std::cerr << "Failed to load font" << std::endl;
  }



  //Load mesh
  if (argc>1) {
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
