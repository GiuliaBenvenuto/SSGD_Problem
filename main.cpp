#include "utilities.h"
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

struct State {
  // Program state
  bool MESH_IS_LOADED;
  // Input
  DrawableTrimesh<> m;  // the input mesh
  uint nverts;          // its #vertices
  vector<vector<uint>> VV;  // its VV relation
  // View
  bool show_m, show_wf;  // show mesh, wire-frame
  // Wireframe
  float wireframe_width;  // Wireframe width
  float wireframe_alpha;  // Wireframe transparency

  State() : MESH_IS_LOADED(false), show_m(true), show_wf(false) {
    // Initialize state with default values
    wireframe_width = 1.0;
    wireframe_alpha = 1.0;
  }
};


//::::::::::::::::::::::::::::::::::::I/O ::::::::::::::::::::::::::::::::::::

void Load_mesh(string filename, GLcanvas & gui, State &gs) {
    gs.m = DrawableTrimesh<>(filename.c_str());
    gs.nverts = gs.m.num_verts();
    gs.VV.resize(gs.nverts); // fill in Vertex-Vertex relation
    for (auto i = 0; i < gs.nverts; i++) gs.VV[i] = gs.m.vert_ordered_verts_link(i);

    // Uncomment the following and adjust parameters if you want a smoother mesh
    // MCF(gs.m, 12, 1e-5, true);

    gs.m.normalize_bbox(); // Rescale mesh to fit [0,1]^3 box
    gs.m.center_bbox();
    gs.m.show_wireframe(gs.show_wf);
    gs.m.show_mesh(gs.show_m);
    gs.m.updateGL();

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

void Setup_GUI_Callbacks(GLcanvas & gui, State &gs) {
  gui.callback_app_controls = [&]() {
    ImGui::SeparatorText("Mesh Operations");

    if (ImGui::Button("Load mesh")) {
      if (gs.MESH_IS_LOADED) {
        ImGui::OpenPopup("Load mesh?");
      } else {
        Load_mesh(gui, gs);
      }
    }

    if (ImGui::BeginPopupModal("Load mesh?", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
      ImGui::Text("This will replace the current mesh.");
      if (ImGui::Button("OK")) {
        Load_mesh(gui, gs);
        ImGui::CloseCurrentPopup();
      }
      ImGui::SameLine();
      if (ImGui::Button("Cancel")) {
        ImGui::CloseCurrentPopup();
      }
      ImGui::EndPopup();
    }

    ImGui::Checkbox("Show Mesh", &gs.show_m);
    if (gs.show_m) {
      gs.m.show_mesh(true);
    } else {
      gs.m.show_mesh(false);
    }

    ImGui::Checkbox("Show Wireframe", &gs.show_wf);
    if (gs.show_wf) {
      gs.m.show_wireframe(true);
    } else {
      gs.m.show_wireframe(false);
    }

    // Wireframe settings
    ImGui::SetNextItemOpen(true, ImGuiCond_Once);
    if (ImGui::TreeNode("Wireframe")) {
      if (ImGui::Checkbox("Show", &gs.show_wf)) gs.m.show_wireframe(gs.show_wf);
      if (ImGui::SliderFloat("Width", &gs.wireframe_width, 1.f, 10.f)) gs.m.show_wireframe_width(gs.wireframe_width);
      // Note: If DrawableTrimesh does not support wireframe transparency directly, this will not work.
      // This line is only useful if you have a method to set wireframe transparency in DrawableTrimesh.
      if (ImGui::SliderFloat("Transparency", &gs.wireframe_alpha, 0.f, 1.f)) {} // Placeholder for transparency control

      ImGui::TreePop();
    }

    if (gs.MESH_IS_LOADED) gs.m.updateGL();
  };
}



//=============================== MAIN =========================================

int main(int argc, char **argv) {

  //SETUP GLOBAL STATE AND GUI:::::::::::::::::::::::
  State gs;
  GLcanvas gui = Init_GUI();
  Setup_GUI_Callbacks(gui,gs);

  //Load mesh
  if (argc>1) {
    string s = "../data/" + string(argv[1]);
    Load_mesh(s,gui,gs);
  }

  // // // GENERATE FIELD:::::::::::::::::::::::::::::::
  // Generate_field(gui,gs);

  // // COMPUTE DISCRETE SCALE SPACE:::::::::::::::::
  // Build_disc_ss(gui,gs);

  // render the mesh
  return gui.launch();
}
