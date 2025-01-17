#include "gc_wrapper.h"

#include <geometrycentral/surface/surface_mesh_factories.h>
#include <cinolib/geometry/vec_mat.h>
// #include <yocto/yocto_cli.h>

using namespace cinolib;
using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace gcHeatWrapper {

// flipout_mesh make_flipout_mesh(
//     const std::vector<vec3i>& triangles, const std::vector<vec3f>& positions) {
//   auto result = flipout_mesh{};

//   // Copy mesh data to compliant format
//   auto polygons           = std::vector<std::vector<size_t>>(triangles.size());
//   auto vertex_coordinates = std::vector<Vector3>(positions.size());
//   for (int i = 0; i < polygons.size(); i++) {
//     auto& tr    = triangles[i];
//     polygons[i] = {(size_t)tr.x, (size_t)tr.y, (size_t)tr.z};
//   }
//   for (int i = 0; i < vertex_coordinates.size(); i++) {
//     auto& p               = positions[i];
//     vertex_coordinates[i] = {p.x, p.y, p.z};
//   }

//   std::tie(result.topology, result.geometry) =
//       makeManifoldSurfaceMeshAndGeometry(polygons, vertex_coordinates);
//   return result;
// }

gc_mesh make_gc_mesh(
    const std::vector<uint>& indices, 
    const std::vector<vec3d>& positions) {
  
  auto result = gc_mesh{};

  // Assuming that every three indices in the 'indices' array form a triangle
  size_t numTriangles = indices.size() / 3;
  auto polygons = std::vector<std::vector<size_t>>(numTriangles);
  auto vertex_coordinates = std::vector<Vector3>(positions.size());

  // Convert the flat index array into the polygon structure expected by makeManifoldSurfaceMeshAndGeometry
  for (size_t i = 0; i < numTriangles; i++) {
    polygons[i] = {indices[3 * i], indices[3 * i + 1], indices[3 * i + 2]};
  }

  // Copy position data, converting from vec3d to Vector3
  for (size_t i = 0; i < positions.size(); i++) {
    auto& p = positions[i];
    vertex_coordinates[i] = {static_cast<float>(p.x()), static_cast<float>(p.y()), static_cast<float>(p.z())};
  }

  // Generate the mesh
  std::tie(result.topology, result.geometry) =
      makeManifoldSurfaceMeshAndGeometry(polygons, vertex_coordinates);

  return result;
}

// Convert vector<uint> to vector<vec3i>
// std::vector<vec3i> convertIndicesToVec3i(const std::vector<uint>& indices) {
//     std::vector<vec3i> triangles;
//     for (size_t i = 0; i < indices.size(); i += 3) {
//         if (i + 2 < indices.size()) {
//             triangles.push_back(vec3i(indices[i], indices[i + 1], indices[i + 2]));
//         }
//     }
//     return triangles;
// }

// // Convert vector<vec3d> to vector<vec3f>
// std::vector<vec3f> convertPositionsToVec3f(const std::vector<vec3d>& positions) {
//     std::vector<vec3f> positions_float;
//     for (const auto& pos : positions) {
//         positions_float.push_back(vec3f(static_cast<float>(pos.x), 
//                                         static_cast<float>(pos.y), 
//                                         static_cast<float>(pos.z)));
//     }
//     return positions_float;
// }
// --------------------------------------------------------------

gc_mesh load_gc_mesh(const std::string& filename) {
  auto result                                = gc_mesh{};
  std::tie(result.topology, result.geometry) = readManifoldSurfaceMesh(
      filename);
  return result;
}

// void shorten_path(
//     FlipEdgeNetwork* edge_network, float angleEPS, bool straightenAtMarked) {
//   // reset counters
//   edge_network->nFlips        = 0;
//   edge_network->nShortenIters = 0;

//   edge_network->EPS_ANGLE                      = angleEPS;
//   edge_network->straightenAroundMarkedVertices = straightenAtMarked;

//   size_t iterLim   = INVALID_IND;
//   double lengthLim = 0.;

//   // START_TIMING(shorten)
//   edge_network->iterativeShorten(iterLim, lengthLim);
//   // FINISH_TIMING_PRINT(shorten)

//   // std::cout << "shortening performed " << edge_network->nShortenIters
//   //           << " iterations, with a total of " << edge_network->nFlips
//   //           << " flips. " << std::endl;
// }

// std::pair<std::unique_ptr<FlipEdgeNetwork>, flipout_path_stats>
// create_path_from_points(ManifoldSurfaceMesh* mesh,
//     VertexPositionGeometry* geometry, int vertex_start, int vertex_end,
//     float angleEPS, bool straightenAtMarked) {
//   if (vertex_start == -1 || vertex_end == -1) return {};

//   auto result =
//       std::pair<std::unique_ptr<FlipEdgeNetwork>, flipout_path_stats>{};

//   // auto  t0           = simple_timer();
//   auto& edge_network = result.first;
//   edge_network       = FlipEdgeNetwork::constructFromDijkstraPath(
//       *mesh, *geometry, mesh->vertex(vertex_start), mesh->vertex(vertex_end));
//   // result.second.initial_guess = elapsed_seconds(t0);

//   // auto t1 = simple_timer();
//   if (edge_network == nullptr) {
//     return result;
//   }
//   edge_network->posGeom = geometry;
//   shorten_path(edge_network.get(), angleEPS, straightenAtMarked);
//   // result.second.shortening = elapsed_seconds(t1);

//   return result;
// }

// std::unique_ptr<FlipEdgeNetwork> make_polyline(ManifoldSurfaceMesh* mesh,
//     VertexPositionGeometry* geometry, const std::vector<int>& vertices,
//     bool closed, bool markInterior) {
//   auto result = std::unique_ptr<FlipEdgeNetwork>(
//       new FlipEdgeNetwork(*mesh, *geometry, {}));
//   result->posGeom = geometry;
//   result->makeDelaunay();

//   auto fullpath = std::vector<Halfedge>{};

//   for (int i = 0; i < vertices.size() - 1; i++) {
//     Vertex vA = result->tri->intrinsicMesh->vertex(vertices[i]);
//     Vertex vB = result->tri->intrinsicMesh->vertex(vertices[i + 1]);
//     if (vA == vB) {
//       throw std::runtime_error("Path starting and ending on the same vertices");
//     }
//     auto path = shortestEdgePath(*result->tri, vA, vB);
//     fullpath.insert(fullpath.end(), path.begin(), path.end());
//   }
//   // Mark control points
//   for (int i = 0; i < vertices.size(); i++) {
//     Vertex v                  = result->tri->intrinsicMesh->vertex(vertices[i]);
//     result->isMarkedVertex[v] = true;
//   }
//   result->addPath(fullpath);

//   return result;
// }

// void subdivide_bezier(FlipEdgeNetwork* control_polygon, int num_subidivisions) {
//   control_polygon->bezierSubdivide(num_subidivisions);
// }

// // Convert to world space positions.
// std::vector<vec3f> path_positions(FlipEdgeNetwork* edge_network) {
//   auto result = std::vector<vec3f>{};

//   // Build a polyline of 3D coordinates
//   std::vector<std::vector<SurfacePoint>> pathPoints =
//       edge_network->getPathPolyline();
//   for (auto& edgePath : pathPoints) {
//     for (const SurfacePoint& p : edgePath) {
//       Vector3 p3d = p.interpolate(edge_network->posGeom->inputVertexPositions);
//       result.push_back({(float)p3d[0], (float)p3d[1], (float)p3d[2]});
//     }
//   }
//   return result;
// }

}  // namespace flipout
