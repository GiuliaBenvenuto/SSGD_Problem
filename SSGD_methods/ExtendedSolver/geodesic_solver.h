#ifndef GEODESIC_SOLVER_H
#define GEODESIC_SOLVER_H

#include <float.h>
#include <iostream>
#include <vector>

using namespace std;

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