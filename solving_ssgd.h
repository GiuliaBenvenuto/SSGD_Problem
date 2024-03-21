#ifndef SOLVING_SSGD_H
#define SOLVING_SSGD_H

#include "SSGD_methods/Graph-based_methods/extended_solver.h"
#include "SSGD_methods/VTP/vtp_wrapper.h"
#include <cinolib/geodesics.h>
#include <cinolib/geometry/vec_mat.h>
#include <cinolib/how_many_seconds.h>
#include <cinolib/meshes/drawable_trimesh.h>
#include <functional>
#include <sys/types.h>

using namespace std;
using namespace cinolib;

// Heat method
ScalarField SSGD_Heat(DrawableTrimesh<> &m,
                      GeodesicsCache &prefactored_matrices,
                      vector<uint> &sources, double &time_heat);

// VTP method
ScalarField SSGD_VTP(DrawableTrimesh<> &m, vector<int> &sources,
                     double &vtp_geodesic_time);

// GeoTangle method
ScalarField SSGD_GeoTangle(DrawableTrimesh<> &m, geodesic_solver &solver,
                           vector<int> &sources,
                           double &geotangle_geodesic_time);

// Edge method
ScalarField SSGD_Edge(DrawableTrimesh<> &m, geodesic_solver &solver,
                      vector<int> &sources, double &edge_geodesic_time);

// Extended method
ScalarField SSGD_Extended(DrawableTrimesh<> &m, geodesic_solver &solver,
                          vector<int> &sources, double &extended_geodesic_time);

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// New Data Structure
////// ABSTRACT CLASS //////
class GeodesicMethod {
public:
  explicit GeodesicMethod() {}
  virtual ~GeodesicMethod() {}

  virtual void load(const std::vector<double> &coords,
                    const std::vector<uint> &tris) = 0;
  virtual void preprocess() = 0;
  virtual void query(const int vid, std::vector<double> &res) = 0;
};

////// INSTANCE OF A SPECIFIC ALGORITHM //////
class ExtendedSolver : public GeodesicMethod {
public:
  ExtendedSolver() {}
  ~ExtendedSolver() {}

  DrawableTrimesh<> m;
  geodesic_solver solver;
  dual_geodesic_solver dual_solver;
  bool dual_solver_computed = false;
  int k = 3;
  void load(const std::vector<double> &coords,
            const std::vector<uint> &tris) override {
    m = DrawableTrimesh(coords, tris);
  }
  void set_k(const int new_k, const bool compute_solver = true) {
    k = new_k;
    if (compute_solver) {
      if (!dual_solver_computed) {
        dual_solver = make_dual_geodesic_solver(m);
        dual_solver_computed = true;
      }
      solver = extended_solver(m, dual_solver, k);
    }
  }

  void preprocess() override {
    dual_solver = make_dual_geodesic_solver(m);
    solver = extended_solver(m, dual_solver, k);
    dual_solver_computed = true;
  }

  void query(const int vid, std::vector<double> &res) override {
    res = compute_geodesic_distances(solver,
                                     {vid}); // TODO change the code so that the
                                             // input is an uint and not a in
  }
}; // TODO: add some checks such as "compute the solver if the query is called
   // but the solver is empty " and stuff like that

#endif