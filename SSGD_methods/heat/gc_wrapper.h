#ifndef FLIPOUT_GC_WRAPPER_H
#define FLIPOUT_GC_WRAPPER_H

#pragma once
#include <geometrycentral/surface/heat_method_distance.h>
#include <geometrycentral/surface/edge_length_geometry.h>
#include <geometrycentral/surface/flip_geodesics.h>
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/mesh_graph_algorithms.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/polygon_soup_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/utilities/timing.h>
#include <cinolib/geometry/vec_mat.h>


//using namespace yocto;
using namespace cinolib;
using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace gcHeatWrapper {

struct gc_mesh {
  std::unique_ptr<ManifoldSurfaceMesh>    topology;
  std::unique_ptr<VertexPositionGeometry> geometry;
};

struct gc_path_stats {
  double initial_guess = 0;
  double shortening    = 0;
};

// vector<vec3i>& triangles -> array of triangles represented by 3 indices of the vertices in the positions array
// vector<vec3f>& positions -> array of 3D positions of the vertices
// I have to change this function such that it takes in input vector<uint> and vector<vec3d>
// flipout_mesh make_flipout_mesh(
//     const std::vector<vec3i>& triangles, const std::vector<vec3f>& positions);
gc_mesh make_gc_mesh(
    const std::vector<uint>& indices, 
    const std::vector<vec3d>& positions);

// Convert vector<uint> to and vector<vec3d> to vec3f
//std::vector<vec3i> convertIndicesToVec3i(const std::vector<uint>& indices);
//std::vector<vec3f> convertPositionsToVec3f(const std::vector<vec3d>& positions);


gc_mesh load_gc_mesh(const std::string& filename);

// Convert to world space positions.
// std::vector<vec3f> path_positions(FlipEdgeNetwork* edge_network);

// std::pair<std::unique_ptr<FlipEdgeNetwork>, flipout_path_stats>
// create_path_from_points(ManifoldSurfaceMesh* mesh,
//     VertexPositionGeometry* geometry, int vertex_start, int vertex_end,
//     float angleEPS = 1e-5, bool straightenAtMarked = true);

// void shorten_path(FlipEdgeNetwork* edge_network, float angleEPS = 1e-5,
//     bool straightenAtMarked = true);

// void subdivide_bezier(FlipEdgeNetwork* control_polygon, int num_subdivisions);

// std::unique_ptr<FlipEdgeNetwork> make_polyline(ManifoldSurfaceMesh* mesh,
//     VertexPositionGeometry* geometry, const std::vector<int>& vertices,
//     bool closed = false, bool markInterior = false);

// inline std::unique_ptr<FlipEdgeNetwork> compute_bezier_curve(
//     const flipout_mesh& mesh, const std::vector<int>& vertices,
//     int subdivisions) {
//   auto bezier = make_polyline(
//       mesh.topology.get(), mesh.geometry.get(), vertices);
//   subdivide_bezier(bezier.get(), subdivisions);
//   return bezier;
// }

}

#endif