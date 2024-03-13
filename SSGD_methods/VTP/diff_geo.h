#ifndef DIFF_GEO_H
#define DIFF_GEO_H
#include "geodesic_algorithm_exact.h"
#include "geodesic_mesh.h"
#include <Eigen/Dense>
#include <cinolib/meshes/drawable_trimesh.h>
using namespace std;
using namespace cinolib;

enum type_of_distance { geodesic, isophotic };

struct patch {
  Eigen::MatrixXd quadric;
  vec3d xu;
  vec3d xv;
  vec3d xuu;
  vec3d xuv;
  vec3d xvv;
  double max_res;
  vector<pair<vec2d, int>> parametric_nbr;
  vector<vec2d> domain_range;
  vector<pair<vec2d, int>> CH;
};

struct geodesic_solver {
  struct graph_edge {
    int node = -1;
    float length = __DBL_MAX__;
    float isophotic_length = __DBL_MAX__;
  };
  vector<vector<graph_edge>> graph = {};
};

struct iVd {
  // map of voronoi verts (triangles)
  vector<int> voronoi_verts = {};
  // vector having the same size of the mesh,
  // vornoi_tags[i]=closest voronoi center of i
  vector<int> voronoi_tags = {};
  // distance field from the voronoi
  // centers
  vector<double> distances = {};

  // regions
  vector<vector<int>> voronoi_regions = {};
  // maximum normal deviation of the vertices
  // in every region w.r.t to the normal at
  // the center
  vector<pair<double, int>> region_normal_deviation = {};

  // PCA of the region
  vector<vector<vec3d>> basis = {};

  // coordinates of the points in the plane
  // defined by PCA
  vector<vector<pair<vec2d, int>>> parametric_nbr = {};

  // maximum in distances
  float R = 0.f;
};

vector<patch> patch_fitting(const DrawableTrimesh<> &m, const int k);

iVd farthest_point_sampling(vector<int> &seeds, const geodesic_solver &solver,
                            const DrawableTrimesh<> &mesh, const int type,
                            const int k);

geodesic_solver compute_geodesic_solver(const DrawableTrimesh<> &m,
                                        const vector<patch> &quadrics);

vector<double> compute_geodesic_distances(const geodesic_solver &solver,
                                          const vector<int> &sources,
                                          const int type);

vector<double> exact_geodesic_distance(const vector<vector<uint>> &triangles,
                                       const vector<vec3d> &positions,
                                       const int &source);

std::tuple<vector<vec3d>, vector<vec3d>>
principal_curvature_field(const vector<patch> &quadrics);

#endif
