#ifndef VTP_WRAPPER
#define VTP_WRAPPER

#include "geodesic_algorithm_exact.h"
#include "geodesic_mesh.h"
#include <cinolib/geometry/vec_mat.h>
using namespace std;
using namespace cinolib;

vector<double> exact_geodesic_distance(const vector<vector<uint>> &triangles,
                                       const vector<vec3d> &positions,
                                       const int &source);
#endif