#include "extended_solver.h"

geodesic_solver extended_solver(const DrawableTrimesh<> &m,
                                const dual_geodesic_solver &solver, int k) {
  geodesic_solver result;
  uint V = m.num_verts();
  result.graph.resize(V);
  for (uint i = 0; i < V; ++i) {
    vector<int> nbr = k_ring(m, i, k);
    result.graph[i].resize(nbr.size());
    for (size_t j = 0; j < nbr.size(); ++j) {
      mesh_point src = get_point_from_vert(m, i);
      mesh_point tgt = get_point_from_vert(m, nbr[j]);
      vector<vec3d> path = shortest_path(src, tgt, m, solver);
      result.graph[i][j].node = nbr[j];
      result.graph[i][j].length = path_length(path);
    }
  }

  return result;
}