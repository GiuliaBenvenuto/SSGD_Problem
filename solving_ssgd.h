#ifndef SOLVING_SSGD_H
#define SOLVING_SSGD_H

#include <cinolib/geometry/vec_mat.h>
#include <cinolib/meshes/drawable_trimesh.h>
#include <cinolib/geodesics.h>
#include "SSGD_methods/VTP/vtp_wrapper.h"

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

#endif