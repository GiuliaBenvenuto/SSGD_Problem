#include "extended_solver.h"

geodesic_solver extended_solver(const DrawableTrimesh<> &m,
                                const dual_geodesic_solver &solver, int k) {
  geodesic_solver result;
  uint V = m.num_verts();
  result.graph.resize(V);

  //cout << "START to computing graph with " << V << " vertices" << endl;
  for (uint i = 0; i < V; ++i) {
    //cout << "Computing vertex " << i << endl;
    vector<int> nbr = k_ring(m, i, k);
    //cout << "Vertex " << i << " has " << nbr.size() << " neighbors" << endl;
    result.graph[i].resize(nbr.size());
    for (size_t j = 0; j < nbr.size(); ++j) {
      auto tic = chrono::high_resolution_clock::now();
      mesh_point src = get_point_from_vert(m, i);
      mesh_point tgt = get_point_from_vert(m, nbr[j]);

      //cout << "Computing SHORTEST PATH between " << i << " and " << nbr[j] << endl;
      vector<vec3d> path = shortest_path(src, tgt, m, solver);
      auto toc = chrono::high_resolution_clock::now();
      // cout time high_resolution_clock in microseconds


      //cout << "SHORTEST PATH computed" << endl;
      result.graph[i][j].node = nbr[j];
      result.graph[i][j].length = path_length(path);
      //cout << "Time: " << chrono::duration_cast<chrono::microseconds>(toc - tic).count() << " microseconds" << endl;
    }
  }
  // // print number of edges
  // uint E = 0;
  // for (uint i = 0; i < V; ++i) {
  //   E += result.graph[i].size();
  // } 
  // cout << "Number of edges: " << E << endl;

  return result;
}


/*
  for (uint i = 0; i < V; ++i) {
    //cout << "Computing vertex " << i << endl;
    vector<int> nbr = k_ring(m, i, k);
    //cout << "Vertex " << i << " has " << nbr.size() << " neighbors" << endl;
    result.graph[i].resize(nbr.size());
    for (size_t j = 0; j < nbr.size(); ++j) {
      if (i > nbr[j]) {
        continue;
      }
      auto tic = chrono::high_resolution_clock::now();
      mesh_point src = get_point_from_vert(m, i);
      mesh_point tgt = get_point_from_vert(m, nbr[j]);

      //cout << "Computing SHORTEST PATH between " << i << " and " << nbr[j] << endl;
      vector<vec3d> path = shortest_path(src, tgt, m, solver);
      auto toc = chrono::high_resolution_clock::now();
      // cout time high_resolution_clock in microseconds


      //cout << "SHORTEST PATH computed" << endl;
      double len = path_length(path);
      result.graph[i][j].node = nbr[j];
      result.graph[i][j].length = len;
      cout << "Time: " << chrono::duration_cast<chrono::microseconds>(toc - tic).count() << " microseconds" << endl;

      if (i < nbr[j]) {
        vector<int> nbr_j = k_ring(m, nbr[j], k);
        if (result.graph[j].size() == 0) {
          result.graph[j].resize(nbr_j.size());
        }
        auto it = find(nbr_j.begin(), nbr_j.end(), i);
        if (it == nbr_j.end()) {
          cout << "Error: vertex " << i << " not found in neighbors of " << nbr[j] << endl;
          exit(1);
        }
        auto entry = distance(nbr_j.begin(), it);
        result.graph[j][entry].node = i;
        result.graph[j][entry].length = len;
      }

      
    }
  }
  */