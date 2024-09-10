#ifndef GEODESIC_SOLVER_H
#define GEODESIC_SOLVER_H
#pragma once
#include <float.h>
#include <iostream>
#include <vector>

#include <cinolib/geometry/vec_mat.h>
#include <cinolib/meshes/drawable_trimesh.h>
using namespace std;
using namespace cinolib;

struct geodesic_solver {
  struct graph_edge {
    int node = -1;
    float length = DBL_MAX;

    // Overload the output operator for graph_edge
    friend ostream &operator<<(ostream &os, const graph_edge &edge) {
      os << "Node: " << edge.node << ", Length: " << edge.length;
      return os;
    }
  };

  vector<vector<graph_edge>> graph = {};

  // Method to print the graph
  void print_graph() {
    for (int i = 0; i < graph.size(); ++i) {
      cout << "Adjacency list of vertex " << i << ": " << endl;
      for (const auto &edge : graph[i]) {
        cout << edge << "; " << endl;
      }
      cout << endl;
    }
  }
};

struct dual_geodesic_solver {
  struct edge {
    int node = -1;
    float length = DBL_MAX;
  };
  vector<array<edge, 3>> graph = {};
};

geodesic_solver make_geodesic_solver(const DrawableTrimesh<> &m,
                                     const bool geo_tangle);

array<vec2d, 3> init_flat_triangle(const DrawableTrimesh<> &m, const uint tid);
array<vec2d, 3> unfold_face(const DrawableTrimesh<> &m, const uint tid,
                            const uint adj, const array<vec2d, 3> &flat_tid);

void update_geodesic_distances(vector<double> &distances,
                               const geodesic_solver &solver,
                               const vector<int> &sources,
                               double max_distance = DBL_MAX);

vector<int> update_distances(vector<double> &distances,
                             const geodesic_solver &solver,
                             const vector<int> &sources,
                             double max_distance = DBL_MAX);

vector<double> compute_geodesic_distances(const geodesic_solver &solver,
                                          const vector<int> &sources);

vector<int> strip_on_dual_graph(const dual_geodesic_solver &solver,
                                const int start, const int end);

vec2d intersect_circles(const vec2d &c2, const double &R2, const vec2d &c1,
                        const double &R1);

inline int vert_offset(const DrawableTrimesh<> &m, const int tid,
                       const int vid) {
  int k = -1;
  vector<uint> vids = m.poly_verts_id(tid);
  for (uint j = 0; j < 3; ++j)
    if (vids[j] == vid)
      return j;

  return k;
};

// ------ Lanthier ------
geodesic_solver compute_fine_graph(DrawableTrimesh<> &m, uint pxedge);
int add_node(geodesic_solver &solver, cinolib::vec3d p);
void add_directed_arc(geodesic_solver &solver, int a, int b, float length);
// add_undirected_arc -> already implemented in connect_nodes()

#endif