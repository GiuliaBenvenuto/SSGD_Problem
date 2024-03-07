#include <cinolib/meshes/meshes.h>
#include <cinolib/geodesics.h>
#include <cinolib/profiler.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>

int main(int argc, char **argv)
{   
    using namespace cinolib;
    const std::string DATA_PATH = "../data/cinolib";
    std::string s = (argc==2) ? std::string(argv[1]) : std::string(DATA_PATH) + "/bunny.obj";
    DrawableTrimesh<> m(s.c_str());
    m.show_wireframe(false);

    GLcanvas gui;
    SurfaceMeshControls<DrawableTrimesh<>> menu(&m, &gui);
    gui.push(&m);
    gui.push(&menu);

    // ---- Position heat sources with CMD + left click ----

    // Create a profiler to measure performance
    Profiler profiler;
    // Create a vector to store the IDs of the heat source vertices on the mesh
    std::vector<uint> sources;
    // Create a GeodesicsCache object to store prefactored matrices
    GeodesicsCache prefactored_matrices;

    // Lambda function that is called when the mouse left button is clicked
    gui.callback_mouse_left_click = [&](int modifiers) -> bool
    {
        // Check if the CMD key (GLFW_MOD_SUPER) is pressed during the mouse click
        if(modifiers & GLFW_MOD_SUPER)
        {
            // Declare a 3D vector to store the unprojected point
            vec3d p;
            // Gets the current cursor position on the screen (a 2D point) and stores it in the "click" variable
            vec2d click = gui.cursor_pos();
            // Try to unproject the 2D cursor position (click) to a 3D point in the scene, storing the result in p
            if(gui.unproject(click, p)) // transform click in a 3d point
            {
                // Call to method pick_vert on an object m (mesh), passing the unprojected point p
                uint vid = m.pick_vert(p);
                // The closest vertex vid to the unprojected point p is identified and added to the sources vector
                sources.push_back(vid);

                // New profiling section named "compute_geodesics"
                profiler.push("compute_geodesics");
                // Compute the geodesic distances from the heat sources to all the other vertices on the mesh
                // The "amortized" part indicates that the function uses some form of optimization to speed up the computation, 
                // likely leveraging the prefactored_matrices which store precomputed data to make repeated geodesic computations more efficient.
                compute_geodesics_amortized(m, prefactored_matrices, sources).copy_to_mesh(m);
                // End the profiling section
                profiler.pop();
                // Show the heat map on the mesh
                m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);
            }
        }
        return false;
    };

    // ---- Reset heat sources with a button ----
    gui.callback_app_controls = [&]()
    {
        if(ImGui::SmallButton("Reset"))
        {
            // Reset heat sources
            sources.clear();
            m.show_vert_color();
        }
    };

    return gui.launch();
}
