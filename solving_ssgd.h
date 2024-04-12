#ifndef SOLVING_SSGD_H
#define SOLVING_SSGD_H

#include "SSGD_methods/Graph-based_methods/extended_solver.h"
#include "SSGD_methods/Trettner/trettner.h"
#include "SSGD_methods/VTP/vtp_wrapper.h"
#include "SSGD_methods/heat/gc_wrapper.h"
#include <geometrycentral/surface/heat_method_distance.h>

#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <cinolib/geodesics.h>
#include <cinolib/geometry/vec_mat.h>
#include <cinolib/how_many_seconds.h>
#include <cinolib/meshes/drawable_trimesh.h>
#include <functional>
#include <sys/types.h>


using namespace std;
using namespace cinolib;
using namespace gcHeatWrapper;

// // Heat method
// ScalarField SSGD_Heat(DrawableTrimesh<> &m,
//                       GeodesicsCache &prefactored_matrices,
//                       vector<uint> &sources, double &time_heat);

// // VTP method
// ScalarField SSGD_VTP(DrawableTrimesh<> &m, vector<int> &sources,
//                      double &vtp_geodesic_time);

// // GeoTangle method
// ScalarField SSGD_GeoTangle(DrawableTrimesh<> &m, geodesic_solver &solver,
//                            vector<int> &sources,
//                            double &geotangle_geodesic_time);

// // Edge method
// ScalarField SSGD_Edge(DrawableTrimesh<> &m, geodesic_solver &solver,
//                       vector<int> &sources, double &edge_geodesic_time);

// // Extended method
// ScalarField SSGD_Extended(DrawableTrimesh<> &m, geodesic_solver &solver,
//                           vector<int> &sources, double &extended_geodesic_time);

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
  virtual void query(const int vid, std::vector<double> &res, ScalarField &sc) = 0;
};

// ---------- VTP ----------
class VTPSolver : public GeodesicMethod {
public:
  VTPSolver() {}
  ~VTPSolver() {}

  DrawableTrimesh<> m;

  void load(const std::vector<double> &coords, const std::vector<uint> &tris) override {
    m = DrawableTrimesh(coords, tris);
  }

  void preprocess() override {}

  void query(const int vid, std::vector<double> &res, ScalarField &sc) override {
    res = exact_geodesic_distance(m.vector_polys(), m.vector_verts(), vid);
    for (auto &value : res) {
      value = 1.0 - value;
    }
    sc = ScalarField(res);
    sc.normalize_in_01();
  }
};


// ---------- Trettner ----------
class TrettnerSolver : public GeodesicMethod {
public:
  TrettnerSolver() {}
  ~TrettnerSolver() {}

  DrawableTrimesh<> m;
  HalfEdge half_edge;
  string mesh_path;
  double time;

  explicit TrettnerSolver(const std::string &path) : mesh_path(path) {}

  void load(const std::vector<double> &coords, const std::vector<uint> &tris) override {
    m = DrawableTrimesh(coords, tris);
  }

  void preprocess() override {}

  void query(const int vid, std::vector<double> &res, ScalarField &sc) override {
    vector<int> vids = {vid};
    half_edge = HEInit(mesh_path, vids);
    sc = distance_field_trettner(half_edge, vids, time); 
  }

};


// ---------- Heat ----------
// class HeatSolver : public GeodesicMethod {
// public:
//   HeatSolver() {}
//   ~HeatSolver() {}

//   DrawableTrimesh<> m;
//   GeodesicsCache prefactored_matrices;
//   bool cache = false;
//   float time_scalar = 1.0;

//   void load(const std::vector<double> &coords, const std::vector<uint> &tris) override {
//     m = DrawableTrimesh(coords, tris);
//   }

//   void preprocess() override {}

//   void set_t(const float new_t) {
//     time_scalar = new_t;
//   }

//   void query(const int vid, std::vector<double> &res, ScalarField &sc) override {
//     if (prefactored_matrices.heat_flow_cache != NULL) {
//       cache = true;
//     }
//     std::vector<uint> vids;
//     vids.push_back(static_cast<uint>(vid));
//     //sc = compute_geodesics_amortized(m, prefactored_matrices, vids);
//     sc = compute_geodesics_amortized(m, prefactored_matrices, vids, COTANGENT, time_scalar);

//     if (cache) {
//     cout << "Heat computation with cache." << endl;
//     } else {
//       cout << "Heat computation without cache." << endl;
//     }
//   }

// };

// ---------- Heat GEOMETRY CENTRAL ----------
class HeatSolver : public GeodesicMethod {
public:
  HeatSolver() {}
  ~HeatSolver() {}

  //DrawableTrimesh<> m;
  //float time_scalar = 1.0;
  gc_mesh gc_m;
  unique_ptr<HeatMethodDistanceSolver> heatSolverGC;
  double time_scalar = 1.0;

  //IntrinsicGeometryInterface intrinsicGeometry;

  void load(const std::vector<double> &coords, const std::vector<uint> &tris) override {
    // from vector<double> to vector<vec3d>
    std::vector<vec3d> converted_coords;
    converted_coords.reserve(coords.size() / 3);

    for (size_t i = 0; i < coords.size(); i += 3) {
        converted_coords.emplace_back(coords[i], coords[i + 1], coords[i + 2]);
    }
    gc_m = make_gc_mesh(tris, converted_coords);

  }

  void preprocess() override {
    // Ensure that the mesh and geometry are loaded and valid
    if (!gc_m.topology || !gc_m.geometry) {
        cerr << "Mesh or geometry not initialized." << endl;
        return;
    }

    // Initialize the heat method distance solver
    heatSolverGC = make_unique<HeatMethodDistanceSolver>(*gc_m.geometry);
  }

  void set_t(const float new_t) {
    time_scalar = new_t;
    // Create the solver considering the new time scalar
    cout << "Time scalar in set_t: " << time_scalar << endl;
    heatSolverGC = make_unique<HeatMethodDistanceSolver>(*gc_m.geometry, time_scalar);
  }

  void query(const int vid, std::vector<double> &res, ScalarField &sc) override {
    // Check if the solver is initialized
    if (!heatSolverGC) {
        cerr << "Heat method solver is not initialized." << endl;
        return;
    }
    // Check if the vertex ID is valid
    if (vid < 0 || vid >= gc_m.topology->nVertices()) {
        cerr << "Invalid vertex ID." << endl;
        return;
    }

    Vertex sourceVertex = Vertex(gc_m.topology.get(), vid);  // Convert int vid to Vertex
    VertexData<double> distances = heatSolverGC->computeDistance(sourceVertex);

    // Process the distances and update the results and ScalarField
    res.clear(); 
    res.reserve(distances.size());
    for (Vertex v : gc_m.topology->vertices()) {
        double value = distances[v];
        value = 1.0 - value;
        res.push_back(value);
    }

    // Update the ScalarField
    sc = ScalarField(res);
    sc.normalize_in_01();
  }

};

// ---------- Geotangle ----------
class GeotangleSolver : public GeodesicMethod {
public:
  GeotangleSolver() {}
  ~GeotangleSolver() {}

  DrawableTrimesh<> m;
  geodesic_solver solver;
  bool solver_computed = false;

  void load(const std::vector<double> &coords, const std::vector<uint> &tris) override {
    m = DrawableTrimesh(coords, tris);
  }

  void preprocess() override {
    solver = make_geodesic_solver(m, true);
    solver_computed = true;
  }

  void query(const int vid, std::vector<double> &res, ScalarField &sc) override {
    res = compute_geodesic_distances(solver, {vid});
    for (auto &value : res) {
      value = 1.0 - value;
    }
    sc = ScalarField(res);
    sc.normalize_in_01();
  }


};

// ---------- Edge ----------
class EdgeSolver : public GeodesicMethod {
public:
  EdgeSolver() {}
  ~EdgeSolver() {}

  DrawableTrimesh<> m;
  geodesic_solver solver;
  bool solver_computed = false;

  void load(const std::vector<double> &coords, const std::vector<uint> &tris) override {
    m = DrawableTrimesh(coords, tris);
  }

  void preprocess() override {
    solver = make_geodesic_solver(m, false);
    solver_computed = true;
  }

  void query(const int vid, std::vector<double> &res, ScalarField &sc) override {

    res = compute_geodesic_distances(solver, {vid});
    for (auto &value : res) {
      value = 1.0 - value;
    }
    sc = ScalarField(res);
    sc.normalize_in_01();
  }
};


// ---------- Extended ----------
class ExtendedSolver : public GeodesicMethod {
public:
  ExtendedSolver() {}
  ~ExtendedSolver() {}

  DrawableTrimesh<> m;
  geodesic_solver solver;
  dual_geodesic_solver dual_solver;
  bool dual_solver_computed = false;
  int k = 3;

  void load(const std::vector<double> &coords, const std::vector<uint> &tris) override {
    m = DrawableTrimesh(coords, tris);
  }
  
  void set_k(const int new_k, const bool compute_solver = true) {
    k = new_k;
    cout << "K in the function: " << k << endl;
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

  void query(const int vid, std::vector<double> &res, ScalarField &sc) override {
    res = compute_geodesic_distances(solver, {vid}); // TODO change the code so that the // input is an uint and not a in
    for (auto &value : res) {
      value = 1.0 - value;
    }
    sc = ScalarField(res);
    sc.normalize_in_01();
  }
}; // TODO: add some checks such as "compute the solver if the query is called
   // but the solver is empty " and stuff like that


// ---------- Lanthier ----------
class LanthierSolver : public GeodesicMethod {
public:
  LanthierSolver() {}
  ~LanthierSolver() {}

  DrawableTrimesh<> m;
  //geodesic_solver solver;
  //dual_geodesic_solver dual_solver;
  //bool dual_solver_computed = false;
  //int k = 3;

  void load(const std::vector<double> &coords, const std::vector<uint> &tris) override {
    m = DrawableTrimesh(coords, tris);
    cout << "Lanthier LOAD" << endl;
  }
  
  void set_k(const int new_k, const bool compute_solver = true) {

  }

  void preprocess() override {
    cout << "Lanthier PREPROCESS" << endl;
  }

  void query(const int vid, std::vector<double> &res, ScalarField &sc) override {

  }

};


#endif