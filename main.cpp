#include "utilities.h"
#include <cinolib/drawable_vector_field.h>
#include <cinolib/drawable_sphere.h>
#include <cinolib/geodesics.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/file_dialog_open.h>
#include <cinolib/gl/file_dialog_save.h>
#include <cinolib/gradient.h>
#include <cinolib/io/write_OBJ.h>
#include <cinolib/mean_curv_flow.h>
#include <fstream>


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

  // SSGD Method
  // enum SSGDMethod { EXACT_POLYHEDRAL, PDE_BASED, GRAPH_BASED } ssgd_method;
  enum SSGDMethod { VTP, HEAT, FMM, GEOTANGLE, EDGE} ssgd_method;

  // Selection of the source vertex
  std::vector<uint> sources;
  GeodesicsCache prefactored_matrices;



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
  }
};



//::::::::::::::::::::::::::::::::::::I/O ::::::::::::::::::::::::::::::::::::
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



//::::::::::::::::::::::::::::::::: GUI ::::::::::::::::::::::::::::::::::::::::::::::::
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
    ImGui::SeparatorText("Files");
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

    // Button for Compute SSGD
    if (ImGui::Button("Compute SSGD")) {
      // Based on the selected SSGD method, perform different actions
      switch (gs.ssgd_method) {

        case State::VTP:
          cout << "Computing SSGD with VTP Method" << endl;
          // Add your code for VTP method here
          break;

        case State::HEAT:
          cout << "Computing SSGD with HEAT Method" << endl;
          // Add your code for HEAT method here
          compute_geodesics_amortized(gs.m, gs.prefactored_matrices, gs.sources).copy_to_mesh(gs.m);
          gs.m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);
          break;

        case State::FMM:
          cout << "Computing SSGD with FMM Method" << endl;
          // Add your code for FMM method here
          break;

        case State::GEOTANGLE:
          cout << "Computing SSGD with GeoTangle Method" << endl;
          // Add your code for GeoTangle method here
          break;

        case State::EDGE:
          cout << "Computing SSGD with Edge Method" << endl;
          // Add your code for Edge method here
          break;

        default:
          cout << "No SSGD method selected" << endl;
          break;
      }
    }
    // Reset button
    if(ImGui::SmallButton("Reset")){
        // Reset heat sources
        gs.sources.clear();
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
        if(modifiers & GLFW_MOD_SUPER) {
            vec3d p;
            vec2d click = gui.cursor_pos();
            if(gui.unproject(click, p)) {
                uint vid = gs.m.pick_vert(p);
                gs.sources.push_back(vid);
                cout << "Source vertex: " << vid << endl;

                // Color the selected vertex in red
                gs.m.vert_data(vid).color = Color::RED();
                gs.m.show_vert_color();

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
  }

  // // // GENERATE FIELD:::::::::::::::::::::::::::::::
  // Generate_field(gui,gs);

  // // COMPUTE DISCRETE SCALE SPACE:::::::::::::::::
  // Build_disc_ss(gui,gs);

  // render the mesh
  return gui.launch();
}
