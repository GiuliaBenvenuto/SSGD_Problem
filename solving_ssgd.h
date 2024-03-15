#ifndef SOLVING_SSGD_H
#define SOLVING_SSGD_H

#include <cinolib/geometry/vec_mat.h>
#include <cinolib/meshes/drawable_trimesh.h>
#include <cinolib/geodesics.h>
#include "SSGD_methods/VTP/vtp_wrapper.h"
#include "SSGD_methods/Graph-based_methods/extended_solver.h"

using namespace std;
using namespace cinolib;

// Heat method
ScalarField SSGD_Heat(DrawableTrimesh<> &m, 
                        GeodesicsCache &prefactored_matrices, 
                        vector<uint> &sources, 
                        double &time_heat);

// VTP method
ScalarField SSGD_VTP(DrawableTrimesh<> &m,
                        vector<int> &sources,
                        double &vtp_geodesic_time);

// GeoTangle method
ScalarField SSGD_GeoTangle(DrawableTrimesh<> &m, 
                        geodesic_solver &solver, 
                        vector<int> &sources,
                        double &geotangle_geodesic_time);

// Edge method
ScalarField SSGD_Edge(DrawableTrimesh<> &m, 
                        geodesic_solver &solver,
                        vector<int> &sources,
                        double &edge_geodesic_time);

// Extended method
ScalarField SSGD_Extended(DrawableTrimesh<> &m,
                        dual_geodesic_solver &solver,
                        vector<int> &sources, 
                        double &extended_geodesic_time);

#endif