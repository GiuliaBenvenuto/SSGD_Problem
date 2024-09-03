#include <cinolib/meshes/drawable_trimesh.h>
#include <cinolib/geometry/vec_mat.h>
#include <cinolib/io/read_write.h>
#include <cinolib/drawable_segment_soup.h>
#include <cinolib/drawable_sphere.h>
#include <cinolib/drawable_vector_field.h>
#include <cinolib/gl/file_dialog_open.h>
#include <cinolib/gl/file_dialog_save.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gradient.h>
#include <cinolib/io/write_OBJ.h>
#include <cinolib/scalar_field.h>
#include <cinolib/color.h>
#include <cinolib/vector_serialization.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

namespace fs = std::filesystem;

struct MeshData {
    std::string filename;
    double avg_degree;
    double edge_density;
    double min_aspect_ratio;
    double min_angle;
    size_t num_verts;
    size_t num_faces;
    double sa_v_ratio;  // Add this line

    MeshData(const std::string& filename, double avg_degree, double edge_density, double min_aspect_ratio, double min_angle, size_t num_verts, size_t num_faces, double sa_v_ratio)
    : filename(filename), avg_degree(avg_degree), edge_density(edge_density), min_aspect_ratio(min_aspect_ratio), min_angle(min_angle), num_verts(num_verts), num_faces(num_faces), sa_v_ratio(sa_v_ratio) {}
};


double triangle_aspect_ratio(const cinolib::vec3d &v0, const cinolib::vec3d &v1, const cinolib::vec3d &v2) {
    double a = (v1 - v0).norm();
    double b = (v2 - v1).norm();
    double c = (v0 - v2).norm();
    double s = (a + b + c) / 2;
    double area = sqrt(s * (s - a) * (s - b) * (s - c));
    double circum_radius = (a * b * c) / (4 * area);
    double in_radius = 2 * area / (a + b + c);
    return circum_radius / in_radius;
}

double triangle_min_angle(const cinolib::vec3d &v0, cinolib::vec3d &v1, cinolib::vec3d &v2) {
    double a = (v1 - v0).norm();
    double b = (v2 - v1).norm();
    double c = (v0 - v2).norm();
    double cosA = (b*b + c*c - a*a) / (2 * b * c);
    double cosB = (a*a + c*c - b*b) / (2 * a * c);
    double cosC = (a*a + b*b - c*c) / (2 * a * b);
    return std::min({acos(cosA), acos(cosB), acos(cosC)}) * 180.0 / M_PI; // convert to degrees
}

double compute_surface_area(const cinolib::DrawableTrimesh<> &mesh) {
    double total_area = 0.0;
    for(size_t pid = 0; pid < mesh.num_polys(); ++pid) {
        std::vector<uint> p_vids = mesh.poly_verts_id(pid);
        const cinolib::vec3d &v0 = mesh.vert(p_vids[0]);
        const cinolib::vec3d &v1 = mesh.vert(p_vids[1]);
        const cinolib::vec3d &v2 = mesh.vert(p_vids[2]);
        double area = 0.5 * ((v1 - v0).cross(v2 - v0)).norm();
        total_area += area;
    }
    return total_area;
}

double compute_volume(const cinolib::DrawableTrimesh<> &mesh) {
    // Assumes the mesh is closed and oriented correctly
    double volume = 0.0;
    for(size_t pid = 0; pid < mesh.num_polys(); ++pid) {
        std::vector<uint> p_vids = mesh.poly_verts_id(pid);
        const cinolib::vec3d &v0 = mesh.vert(p_vids[0]);
        const cinolib::vec3d &v1 = mesh.vert(p_vids[1]);
        const cinolib::vec3d &v2 = mesh.vert(p_vids[2]);
        volume += (v0.dot(v1.cross(v2))) / 6.0;
    }
    return std::fabs(volume);
}

MeshData analyze_mesh(const std::string &filename) {
    cinolib::DrawableTrimesh<> mesh(filename.c_str()); // Load mesh

    double avg_degree = 0;
    for (size_t vid = 0; vid < mesh.num_verts(); ++vid) {
        avg_degree += mesh.vert_valence(vid);
    }
    avg_degree /= mesh.num_verts();

    double edge_density = static_cast<double>(mesh.num_edges()) / mesh.num_verts();
    double min_aspect_ratio = std::numeric_limits<double>::max();
    double min_angle = 180.0;

    for (size_t pid = 0; pid < mesh.num_polys(); ++pid) {
        std::vector<uint> p_vids = mesh.poly_verts_id(pid);
        cinolib::vec3d v0 = mesh.vert(p_vids[0]);
        cinolib::vec3d v1 = mesh.vert(p_vids[1]);
        cinolib::vec3d v2 = mesh.vert(p_vids[2]);
        double aspect_ratio = triangle_aspect_ratio(v0, v1, v2);
        double angle = triangle_min_angle(v0, v1, v2);
        min_aspect_ratio = std::min(min_aspect_ratio, aspect_ratio);
        min_angle = std::min(min_angle, angle);
    }

    double surface_area = compute_surface_area(mesh);
    double volume = compute_volume(mesh);
    double sa_v_ratio = surface_area / volume;

    return MeshData(filename, avg_degree, edge_density, min_aspect_ratio, min_angle, mesh.num_verts(), mesh.num_polys(), sa_v_ratio);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <folder_path>" << std::endl;
        return 1;
    }

    std::string folder_path = argv[1];
    std::vector<MeshData> mesh_data_list;

    for (const auto &entry : fs::directory_iterator(folder_path)) {
        if (entry.path().extension() == ".obj") {
            std::string mesh_path = entry.path().string();
            std::cout << "Analyzing mesh: " << mesh_path << std::endl;
            mesh_data_list.push_back(analyze_mesh(mesh_path));
        }
    }

    // Sort the mesh data by the number of vertices
    std::sort(mesh_data_list.begin(), mesh_data_list.end(), [](const MeshData& a, const MeshData& b) {
        return a.num_verts < b.num_verts;
    });

    // Write sorted data to file
    std::ofstream out("all_meshes_metrics.txt");
    if (!out) {
        std::cerr << "Failed to create output file." << std::endl;
        return 1;
    }

    for (const auto& data : mesh_data_list) {
        out << "Mesh: " << data.filename << std::endl;
        out << "Number of Vertices: " << data.num_verts << std::endl;
        out << "Number of Faces: " << data.num_faces << std::endl;
        out << "Average Degree: " << data.avg_degree << std::endl;
        out << "Edge Density: " << data.edge_density << std::endl;
        out << "Minimum Aspect Ratio: " << data.min_aspect_ratio << std::endl;
        out << "Minimum Angle: " << data.min_angle << std::endl;
        out << "Surface Area to Volume Ratio: " << data.sa_v_ratio << std::endl;  // Add this line
        out << std::endl;  // Empty line for readability
    }

    out.close();
    return 0;
}
