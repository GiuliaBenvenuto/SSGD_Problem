#include "vtp_wrapper.h"

std::pair<vector<double>, vector<unsigned>>
exact_geodesic_wrapper(const vector<vector<uint>> &triangles,
                       const vector<vec3d> &positions) {
  int V = (int)positions.size();
  int F = (int)triangles.size();
  vector<double> points(3 * V);
  vector<uint> faces(3 * F);
  vector<double> f(V);

  for (int i = 0; i < V; ++i) {
    points[3 * i] = positions[i].x();
    points[3 * i + 1] = positions[i].y();
    points[3 * i + 2] = positions[i].z();
  }
  for (int i = 0; i < F; ++i) {
    faces[3 * i] = triangles[i][0];
    faces[3 * i + 1] = triangles[i][1];
    faces[3 * i + 2] = triangles[i][2];
  }

  return std::make_pair(points, faces);
}
vector<double> exact_geodesic_distance(const vector<vector<uint>> &triangles,
                                       const vector<vec3d> &positions,
                                       const int &source) {
  int V = (int)positions.size();
  int F = (int)triangles.size();
  vector<double> f;
  f.resize(V);
  auto [points, faces] = exact_geodesic_wrapper(triangles, positions);
  geodesic_VTP::Mesh mesh;
  mesh.initialize_mesh_data(points, faces);
  geodesic_VTP::GeodesicAlgorithmExact algorithm(&mesh);
  algorithm.propagate(source);
  vector<geodesic_VTP::Vertex> verts = mesh.vertices();
  for (int j = 0; j < V; ++j) {
    geodesic_VTP::Vertex v = verts[j];
    auto value = v.geodesic_distance();
    f[j] = value;
  }

  return f;
}