#ifndef TRET
#define TRET

#include <cinolib/geometry/vec_mat.h>
#include <cinolib/meshes/drawable_trimesh.h>
#include <cinolib/geodesics.h>
#include <vector>
#include <cmath>
using namespace std;
using namespace cinolib;

inline float pow2(float f) { return f * f; }
inline float sqrt_sat(float f) { return f <= 0.0f ? 0.0f : std::sqrt(f); }
inline float fast_abs(float f) { return f < 0 ? -f : f; }

struct triangle
{
    int v0 = 0;
    int v1 = 0;
    int v2 = 0;
};

struct pos3
{
    float x = 0;
    float y = 0;
    float z = 0;

    pos3 operator+(pos3 const& r) const { return {x + r.x, y + r.y, z + r.z}; }
    pos3 operator/(float f) const { return {x / f, y / f, z / f}; }
};

inline float dist(pos3 a, pos3 b)
{
    auto dx = a.x - b.x;
    auto dy = a.y - b.y;
    auto dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

struct HalfEdge {
    vector<pos3> vertex_pos;  // Use vec3d from cinolib for vertices
    vector<triangle> faces;  // Ensure 'triangle' is defined or included
    vector<int> halfedge_to_vertex;
    vector<int> halfedge_to_face;
    vector<int> halfedge_to_next;
    vector<int> halfedge_to_prev;
    vector<int> vertex_to_outgoing_halfedge;
    vector<float> edge_lengths;
    float avg_edge_length;
};


// Declare the HEInit function
HalfEdge HEInit(const string &file, 
                vector<int> &sources);

// Declare the compute_distance_field function
ScalarField distance_field_trettner(const HalfEdge &mesh, 
                                    const vector<int> &sources,
                                    double &trettner_geodesic_time);


#endif