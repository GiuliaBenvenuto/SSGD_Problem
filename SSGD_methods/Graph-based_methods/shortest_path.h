#ifndef SHORTEST_PATH_H
#define SHORTEST_PATH_H

#include "geodesic_solver.h"
#include <cinolib/geometry/vec_mat.h>
#include <cinolib/meshes/drawable_trimesh.h>
using namespace std;
using namespace cinolib;

struct mesh_point {
  int tid = -1;
  vec3d bary = vec3d(0, 0, 0);
  mesh_point() {
    tid = -1;
    bary = vec3d(0., 0., 0);
  }
  mesh_point(const int id, const vec3d b) {
    tid = id;
    bary = b;
  }
};

inline std::ostream &operator<<(std::ostream &os, const mesh_point &p) {
  os << "mesh_point(" << p.tid << ","
     << "vec3d(" << p.bary.x() << "," << p.bary.y() << "," << p.bary.z()
     << "))";
  return os;
}

struct funnel_point {
  int face = 0;
  vec2d pos;
};

vector<int> k_ring(const DrawableTrimesh<> &m, const int vid, const int k,
                   const bool boundary_only = true);

mesh_point get_point_from_vert(const DrawableTrimesh<> &m, const int vid);

dual_geodesic_solver make_dual_geodesic_solver(const DrawableTrimesh<> &m);

vector<vec3d> shortest_path(mesh_point &src, mesh_point &tgt,
                            const DrawableTrimesh<> &m,
                            const dual_geodesic_solver &solver);

double shortest_path_polyhedral_distance(mesh_point &src, mesh_point &tgt,
                                         const DrawableTrimesh<> &m,
                                         const dual_geodesic_solver &solver);

double path_length(const vector<vec3d> &path);

template <typename T>
inline T interpolate(const T &p, const T &q, const double alpha) {
  return alpha * p + (1 - alpha) * q;
}
template <typename T>
inline T interpolate(const T &p0, const T &p1, const T &p2,
                     const vec3d &alpha) {
  return alpha[0] * p0 + alpha[1] * p1 + alpha[2] * p2;
}
template <typename T>
inline T interpolate(const array<T, 3> &v, const vec3d &alpha) {
  return alpha[0] * v[0] + alpha[1] * v[1] + alpha[2] * v[2];
}
template <typename T>
inline T interpolate(const vector<T> &v, const vec3d &alpha) {
  return alpha[0] * v[0] + alpha[1] * v[1] + alpha[2] * v[2];
}

#endif
