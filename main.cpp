// #include "utilities.h"
#include <cinolib/drawable_vector_field.h>
#include <cinolib/drawable_segment_soup.h>
#include <cinolib/drawable_sphere.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/file_dialog_open.h>
#include <cinolib/gl/file_dialog_save.h>
#include <cinolib/gradient.h>
#include <cinolib/meshes/drawable_tetmesh.h>
#include <cinolib/vector_serialization.h>
#include <cinolib/io/write_OBJ.h>
#include <cinolib/mean_curv_flow.h>
#include <Eigen/SparseCholesky>
#include <fstream>

// SSDG with Heat method
#include <cinolib/geodesics.h>

// SSDG with VTP method
#include "SSGD_methods/VTP/diff_geo.h"
#include "SSGD_methods/VTP/diff_geo.cpp"

// SSGD with GeoTangle 
//#include "SSGD_methods/GeoTangle/GeoTangle.cpp"
#include "SSGD_methods/GeoTangle/GeoTangle.cpp"

// SSGD with Edge method
#include "SSGD_methods/Edge/edge.cpp"

using namespace std;
using namespace cinolib;

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
  // vec3d vec_color;
  Color vec_color;

  Color vert_color;
  Color poly_color;

  // SSGD Method
  // enum SSGDMethod { EXACT_POLYHEDRAL, PDE_BASED, GRAPH_BASED } ssgd_method;
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
    vecfield_size = 2;
    // vec_color = vec3d(1.0, 0.0, 0.0);
    vec_color = Color(1.0, 0.0, 0.0, 1.0);

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
  // Method inside: #include <cinolib/geodesics.h>
  compute_geodesics_amortized(m, prefactored_matrices, sources).copy_to_mesh(m);
  m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);
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
  gui.side_bar_width = 0.25;
  gui.show_side_bar = true;
  return gui;
}

void Setup_GUI_Callbacks(GLcanvas & gui, State &gs)
{
  gui.callback_app_controls = [&]() {
    // Files
    ImGui::SeparatorText("Single Source Geodesic Distance Computation");
    ImGui::SetNextItemOpen(true, ImGuiCond_Once);
    if (ImGui::TreeNode("IO")) {
      if (ImGui::Button("Load mesh")) {
        if (gs.MESH_IS_LOADED) {
          ImGui::OpenPopup("Load mesh?");
        } else {
          Load_mesh(gui,gs);
        }
      }
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
    }
    

    // Wireframe settings
    ImGui::SetNextItemOpen(true, ImGuiCond_Once);
    if (ImGui::TreeNode("Wireframe")) {
      if (ImGui::Checkbox("Show", &gs.SHOW_WIREFRAME)) gs.m.show_wireframe(gs.SHOW_WIREFRAME);
      if (ImGui::SliderFloat("Width", &gs.wireframe_width, 1.f, 10.f)) gs.m.show_wireframe_width(gs.wireframe_width);
      if (ImGui::SliderFloat("Transparency", &gs.wireframe_alpha, 0.f, 1.f)) gs.m.show_wireframe_transparency(gs.wireframe_alpha);
      ImGui::TreePop();
    }

    // Mesh shading
    ImGui::SetNextItemOpen(true, ImGuiCond_Once);
    if (ImGui::TreeNode("Shading")) {
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
      ImGui::TreePop();
    }

    ImGui::SetNextItemOpen(true, ImGuiCond_Once);
    if (ImGui::TreeNode("Colors")) {
        if (ImGui::BeginTable("Color by:", 2)) {
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            if (ImGui::RadioButton("Vert", gs.m.drawlist.draw_mode & DRAW_TRI_VERTCOLOR)) {
                gs.m.show_vert_color();
            }
            ImGui::TableNextColumn();
            if (ImGui::ColorEdit4("Vertex Color", gs.vert_color.rgba)) {
                // Apply vertex color
                /* for (uint vid = 0; vid < gs.m.num_verts(); ++vid) {
                    gs.m.vert_data(vid).color = gs.vert_color;
                } */
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
                // Apply polygon color
                /* for (uint fid = 0; fid < gs.m.num_polys(); ++fid) {
                    gs.m.poly_data(fid).color = gs.poly_color;
                } */
                gs.m.poly_set_color(gs.poly_color);
                gs.m.show_poly_color();
                gs.m.updateGL();
            }

            ImGui::EndTable();
        }
        ImGui::TreePop();
    }




    // SSGD Method
    ImGui::SetNextItemOpen(true, ImGuiCond_Once);
    if (ImGui::TreeNode("SSGD Method")) {
      // Exact Polyhedral Methods
      ImGui::SeparatorText("Exact Polyhedral Methods");
      if(ImGui::RadioButton("VTP ", gs.ssgd_method == State::VTP)) {
        gs.ssgd_method = State::VTP;
        cout << "VTP Method" << endl;
      }
      // PDE-Based Methods
      ImGui::SeparatorText("PDE-Based Methods");
      if(ImGui::RadioButton("Heat  ", gs.ssgd_method == State::HEAT)) {
        gs.ssgd_method = State::HEAT;
        cout << "HEAT Method" << endl;
      }
      if(ImGui::RadioButton("FMM  ", gs.ssgd_method == State::FMM)) {
        gs.ssgd_method = State::FMM;
        cout << "FMM Method" << endl;
      }
      // Graph-Based Methods
      ImGui::SeparatorText("Graph-Based Methods");
      if(ImGui::RadioButton("GeoTangle  ", gs.ssgd_method == State::GEOTANGLE)) {
        gs.ssgd_method = State::GEOTANGLE;
        cout << "GeoTangle Method" << endl;
      }
      if(ImGui::RadioButton("Edge  ", gs.ssgd_method == State::EDGE)) {
        gs.ssgd_method = State::EDGE;
        cout << "Edge Method" << endl;
      }

      ImGui::TreePop();
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
            cout << "No source selected" << endl;
            ImGui::OpenPopup("Warning");
        } else {

          switch (gs.ssgd_method) {
      
            case State::VTP: {
              cout << "Computing SSGD with VTP Method" << endl;
              SSGD_VTP(gs.m, gs.solver, gs.field_data, gs.field, gs.voronoi_centers);
              break;
            }

            case State::HEAT: {
              cout << "Computing SSGD with HEAT Method" << endl;
              SSGD_Heat(gs.m, gs.prefactored_matrices, gs.sources);
              break;
            }

            case State::FMM: {
              cout << "Computing SSGD with FMM Method" << endl;
              // Add your code for FMM method here
              break;
            }

            case State::GEOTANGLE: {
              cout << "Computing SSGD with GeoTangle Method" << endl;
              SSGD_GeoTangle(gs.m, gs.solver_geo, gs.field_geo, gs.sources_geo);
              break;
            }

            case State::EDGE: {
              cout << "Computing SSGD with Edge Method" << endl;
              SSGD_Edge(gs.m, gs.solver_edge, gs.field_edge, gs.sources_edge);
              break;
            }

            default:
              cout << "No SSGD method selected" << endl;
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

    /* Vector field
    ImGui::SetNextItemOpen(true, ImGuiCond_Once);
    if (ImGui::TreeNode("Vector Field")) {
      if (ImGui::Checkbox("Show Vector Field", &gs.show_vecfield)) {
        cout << "Show vector field: " << gs.show_vecfield << endl;
        if (gs.show_vecfield) {
            cout << "Qui: " << endl;
            if (gs.vec_field.size() == 0) { // Initialize the vector field if it's not already
                cout << "Vector field size before initialization: " << gs.vec_field.size() << endl;
                cout << "Vector field arrow size: " << gs.vecfield_size << endl;

                gs.vec_field = DrawableVectorField(gs.m);
                ScalarField f(gs.m.serialize_uvw(3)); // Assuming U_param is defined somewhere relevant
  
                cout << "Scalar field size: " << f.size() << endl;

                gs.vec_field = gradient_matrix(gs.m) * f;
                gs.vec_field.normalize();

                cout << "Prova" << float(gs.m.edge_avg_length()) * gs.vecfield_size << endl;
                gs.vec_field.set_arrow_size(float(gs.m.edge_avg_length()) * gs.vecfield_size);
                //gs.vec_field.set_arrow_color(Color(gs.vec_color.x(), gs.vec_color.y(), gs.vec_color.z()));
                gs.vec_field.set_arrow_color(gs.vec_color);
                cout << "color" << gs.vec_color << endl;

                cout << "Vector field size after initialization: " << gs.vec_field.size() << endl;

            }
            cout << "Push" << endl;
            gui.push(&gs.vec_field, false);
        } else {
            cout << "Pop" << endl;
            gui.pop(&gs.vec_field);
        }
      }
      // Slider to adjust the vector field arrow size
      if (ImGui::SliderFloat("Size", &gs.vecfield_size, 0.1f, 5.f)) {
          gs.vec_field.set_arrow_size(float(gs.m.edge_avg_length()) * gs.vecfield_size *5);
          cout << "Arrow size: " << float(gs.m.edge_avg_length()) * gs.vecfield_size << endl;
      }

      if (ImGui::ColorEdit4("Color##vec", (float*)&gs.vec_color)) {
          gs.vec_field.set_arrow_color(gs.vec_color);
      }

      ImGui::TreePop();
    }*/
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
                std::cout << "Selected vid VTP = " << selected_vid << std::endl;

                // GeoTangle sources
                gs.sources_geo.push_back(selected_vid);
                std::cout << "Selected vid GEO = " << selected_vid << std::endl;

                // Edge sources
                gs.sources_edge.push_back(selected_vid);
                std::cout << "Selected vid EDGE = " << selected_vid << std::endl;


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

  //Load mesh
  if (argc>1) {
    string s = "../data/" + string(argv[1]);
    Load_mesh(s, gui, gs);
  } else {
    string s = "../data/cinolib/3holes.obj";
    Load_mesh(s, gui, gs);
  }

  // // // GENERATE FIELD:::::::::::::::::::::::::::::::
  // Generate_field(gui,gs);

  // // COMPUTE DISCRETE SCALE SPACE:::::::::::::::::
  // Build_disc_ss(gui,gs);

  // render the mesh
  return gui.launch();
}
