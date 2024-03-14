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
  };
  vector<vector<graph_edge>> graph = {};
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

#endif