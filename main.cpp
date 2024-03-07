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
  // Program state ::::::::::::::::::::::::::::::::::::::::::::::::::::::
  bool MESH_IS_LOADED, FIELD_IS_PRESENT, SCALE_SPACE_IS_PRESENT, EIGENFUNCTIONS_COMPUTED;
  // Input
  DrawableTrimesh<> m;            // the input mesh
  uint nverts;                    // its #vertices
  vector<vector<uint>> VV;        // its VV relation
  // Field
  Eigen::VectorXd f;              // base field on m
  float f_clamp[2];               // clamp limits of base field
  vector<double> eigenfunctions;  // eigenfunctions of f as computed by SPECTRA
  // Scale-space
  vector<vector<double>> fields;  // the discrete scale-space
  vector<vector<int>> critical;   // critical points at all levels
  vector<DrawableSphere> points;  // spheres for rendering critical points
  float point_size;               // base radius of spheres
  // GUI state ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // Field
  int current_field_method;       // current method to generate the field
  bool normalize_f;               // base field is normalized
  int max_eigenfunctions, selected_eigenfunction; // parameters for eigenfunction generator
  // Scale space
  bool normalize;                 // fields are normalized at all levels during diffusion
  int nlevels;                    // # levels in the scale-space 
  int method;                     // 1 diffusion; 2 smoothness optiization
  int diff_progression;           // 1 linear - 2 exponential
  float diff_lambda;              // base time step of diffusion
  int diff_stride;                // stride for linear method
  float diff_mult;                // multiplicator for exponential method
  int opt_progression;            // 1 linear - 2 exponential
  float opt_w;                    // base regularization of optimization
  int opt_stride;                 // stride for linear method
  float opt_div;                  // divisor for exponential method
  // View
  bool show_m, show_wf, show_cp;  // show mesh, wire-frame, critical points
  int show_field;                 // field to show: 0 no field; 1 base field; 2 scale-space field
  float point_multiplier;         // multiplier of sphere radius
  vector<float*> clamp_limits;    // clamp limits for field visualization at all levels
  float curr_clamp[2];            // current clamp limits
  int selected_entry;             // selected level for visualization

  State() {
    MESH_IS_LOADED = FIELD_IS_PRESENT = SCALE_SPACE_IS_PRESENT = false;
    EIGENFUNCTIONS_COMPUTED = false;
    // field generation
    normalize_f = false;
    max_eigenfunctions = 100;
    selected_eigenfunction = 1;
    f_clamp[0] = 0; f_clamp[1] = 1;
    // scale-space
    normalize = true;
    nlevels = 300;
    method = 1;           //diffusion
    diff_progression = 1;  //linear
    diff_lambda = 0.0001;
    diff_stride = 10;
    diff_mult = 1.05;
    opt_progression = 1;  //linear
     opt_w = 1e6;
    opt_stride = 10;
    opt_div = 1.05;
    // view
    show_m = true;
    show_cp = show_wf = false;
    show_field = 0;
    point_multiplier = 1.0;
    f.resize(2); f(0)=0.0; f(1)=1.0;  // init f to support clamp limits in gui
    curr_clamp[0]=0.0; curr_clamp[1]=1.0;
    selected_entry = 0;
  }
};


//::::::::::::::::::::::::::::::::::::GUI utilities ::::::::::::::::::::::::::::::::::::
// functions to  render vertices as spheres
inline void remove_points(const vector<DrawableSphere> &cp, GLcanvas &gui) {
  for (auto &point : cp) gui.pop(&point);
}


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
  gs.m.show_wireframe(gs.show_wf);
  gs.m.show_mesh(gs.show_m);    
  gs.m.updateGL();  
  gs.point_size = gs.m.edge_avg_length()/2; // set initial radius of spheres for critical points
  if (!gs.MESH_IS_LOADED) {
    gui.push(&gs.m);
    gs.MESH_IS_LOADED = true;
  }

  // reset field and scale-space
  gs.f.resize(2); gs.f(0)=0.0; gs.f(1)=1.0; // init f to support clamp limits in gui
  gs.curr_clamp[0]=0.0; gs.curr_clamp[1]=1.0;
  gs.show_field = 0;

  if (gs.show_cp) {
    remove_points(gs.points,gui);
    gs.show_cp = false;
  }

  gs.clamp_limits = vector<float*>(1); // only entry zero for the base field
  gs.clamp_limits[0] = new float[2];

  gs.FIELD_IS_PRESENT = gs.SCALE_SPACE_IS_PRESENT = gs.EIGENFUNCTIONS_COMPUTED = false; 
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
    ImGui::SeparatorText("Files");

    ImGui::SeparatorText("Field");
    
    ImGui::SeparatorText("Scale-space");
    
    ImGui::SeparatorText("View");
  };

  gui.callback_mouse_left_click = [&](int modifiers) -> bool {
    if (modifiers & GLFW_MOD_SHIFT) {
      vec3d p;
      vec2d click = gui.cursor_pos();
      if (gui.unproject(click, p)) {
        uint vid = gs.m.pick_vert(p);
        cout << "Picked vertex " << vid << " field value " 
              << std::setprecision(20) << gs.fields[gs.selected_entry][vid] 
              << " at level " << gs.selected_entry << endl;
        // m.vert_data(vid).color = Color::RED();
        // m.updateGL();
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
