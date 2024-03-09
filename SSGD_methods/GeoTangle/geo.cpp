// #include "diff_geo.h"
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <algorithm>
#include <cinolib/cino_inline.h>
#include <cinolib/geometry/triangle_utils.h>
#include <cinolib/vertex_mass.h>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <vector>
using namespace cinolib;
using namespace Eigen;
using namespace std;
using namespace std::complex_literals;

#include "../../libs/ann_1.1.2/include/ANN/ANN.h"
#include "../../libs/ann_1.1.2/include/ANN/ANNx.h"
#include "../../libs/ann_1.1.2/include/ANN/ANNperf.h"


// Useful: ---------------------------- Graph Construction ----------------------------
void connect_nodes(geodesic_solver &solver, int a, int b, float length) {
  solver.graph[a].push_back({b, length});
  solver.graph[b].push_back({a, length});
}

double opposite_nodes_arc_length(const vector<vec3d> &positions, int a, int c,
                                 const vec2i &edge) {
  // Triangles (a, b, d) and (b, d, c) are connected by (b, d) edge
  // Nodes a and c must be connected.

  auto b = edge.x(), d = edge.y();
  auto ba = positions[a] - positions[b];
  auto bc = positions[c] - positions[b];
  auto bd = positions[d] - positions[b];
  ba.normalize();
  bd.normalize();
  auto cos_alpha = ba.dot(bd);
  auto cos_beta = bc.dot(bd);
  auto sin_alpha = sqrt(max(0.0, 1 - cos_alpha * cos_alpha));
  auto sin_beta = sqrt(max(0.0, 1 - cos_beta * cos_beta));

  // cos(alpha + beta)
  auto cos_alpha_beta = cos_alpha * cos_beta - sin_alpha * sin_beta;
  if (cos_alpha_beta <= -1)
    return DBL_MAX;

  // law of cosines (generalized Pythagorean theorem)
  ba = positions[a] - positions[b];
  bc = positions[c] - positions[b];
  bd = positions[d] - positions[b];
  auto len =
      ba.dot(ba) + bc.dot(bc) - ba.norm() * bc.norm() * 2 * cos_alpha_beta;

  if (len <= 0)
    return DBL_MAX;
  else
    return sqrt(len);
}

static void connect_opposite_nodes(geodesic_solver &solver,
                                   const vector<vec3d> &positions,
                                   const vector<uint> &tr0,
                                   const vector<uint> &tr1, const vec2i &edge) {
  auto opposite_vertex = [](const vector<uint> &tr, const vec2i &edge) -> int {
    for (auto i = 0; i < 3; ++i) {
      if (tr[i] != edge.x() && tr[i] != edge.y())
        return tr[i];
    }
    return -1;
  };

  auto v0 = opposite_vertex(tr0, edge);
  auto v1 = opposite_vertex(tr1, edge);
  if (v0 == -1 || v1 == -1)
    return;
  auto length = opposite_nodes_arc_length(positions, v0, v1, edge);
  connect_nodes(solver, v0, v1, length);
}

geodesic_solver make_geodesic_solver(const DrawableTrimesh<> &m) {
  auto solver = geodesic_solver{};
  solver.graph.resize(m.num_verts());
  for (auto face = 0; face < m.num_polys(); face++) {
    for (auto k = 0; k < 3; k++) {
      auto a = m.poly_vert_id(face, k);
      auto b = m.poly_vert_id(face, (k + 1) % 3);

      // connect mesh edges
      auto len = (m.vert(a) - m.vert(b)).norm();
      if (a < b)
        connect_nodes(solver, a, b, len);

      // connect opposite nodes
      auto neighbor = m.adj_p2p(face)[k];
      if (face < neighbor) {
        connect_opposite_nodes(solver, m.vector_verts(), m.adj_p2v(face),
                               m.adj_p2v(neighbor), vec2i{(int)a, (int)b});
      }
    }
  }
  return solver;
}
// ---------------------------------------------------
