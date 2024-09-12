#ifndef SOLVING_SSGD_H
#define SOLVING_SSGD_H

#pragma once
#include "SSGD_methods/Graph-based_methods/extended_solver.h"
#include "SSGD_methods/Trettner/trettner.h"
#include "SSGD_methods/VTP/vtp_wrapper.h"
#include "SSGD_methods/heat/gc_wrapper.h"
#include <cinolib/geodesics.h>
#include <cinolib/geometry/vec_mat.h>
#include <cinolib/how_many_seconds.h>
#include <cinolib/meshes/drawable_trimesh.h>
#include <functional>
#include <geometrycentral/surface/fast_marching_method.h>
#include <geometrycentral/surface/heat_method_distance.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <sys/types.h>

// Matlab
#include <MatlabDataArray.hpp>
#include <MatlabEngine.hpp>

using namespace std;
using namespace cinolib;
using namespace gcHeatWrapper;

using namespace matlab::engine;
using namespace matlab::data;

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// New Data Structure
////// ABSTRACT CLASS //////

// ------ Helper functions ------
inline vector<double> extract_coords(const DrawableTrimesh<> &mesh) {
  vector<double> coords;
  auto verts = mesh.vector_verts();
  for (const auto &vert : verts) {
    coords.push_back(vert.x());
    coords.push_back(vert.y());
    coords.push_back(vert.z());
  }
  return coords;
}

inline vector<uint> extract_tris(const DrawableTrimesh<> &mesh) {
  vector<uint> tris;
  auto polys = mesh.vector_polys();
  for (const auto &poly : polys) {
    for (auto vid : poly) {
      tris.push_back(vid);
    }
  }
  return tris;
}

inline bool check_sp_polyhedral_distance(const DrawableTrimesh<> &m,
                                         const geodesic_solver &solver) {
  for (size_t i = 0; i < solver.graph.size(); ++i) {
    vector<double> d =
        exact_geodesic_distance(m.vector_polys(), m.vector_verts(), i);
    for (size_t j = 0; j < solver.graph[i].size(); ++j) {
      double len = solver.graph[i][j].length;
      int nei = solver.graph[i][j].node;
      double real_len = d[nei];
      if (std::abs(len - real_len) > 0.1) {
        assert(false);
        return false;
      }
    }
  }
  return true;
}

class GeodesicMethod {
public:
  explicit GeodesicMethod() {}
  virtual ~GeodesicMethod() {}

  virtual void load(DrawableTrimesh<> *mesh) = 0;
  virtual void preprocess() = 0;
  virtual void query(const int vid, std::vector<double> &res) = 0;
};

// ---------- VTP ----------
class VTPSolver : public GeodesicMethod {
public:
  VTPSolver() {}
  ~VTPSolver() {}

  DrawableTrimesh<> *m;

  void load(DrawableTrimesh<> *mesh) override { m = mesh; }

  void preprocess() override {}

  void query(const int vid, std::vector<double> &res) override {
    // cout << "VTP query started... " << endl;
    res.clear(); // Clear previous results

    auto start = std::chrono::high_resolution_clock::now();
    res = exact_geodesic_distance(m->vector_polys(), m->vector_verts(), vid);

    // for (int i = 0; i < res.size(); i++) {
    //   cout << i << "," << res[i] << endl;
    // }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    // cout << "VTP computation: " << elapsed.count() << " s" << endl;
  }
};

// ---------- Trettner ----------
class TrettnerSolver : public GeodesicMethod {
public:
  TrettnerSolver() {}
  ~TrettnerSolver() {}

  HalfEdge half_edge;
  string mesh_path;
  double time;

  void load(DrawableTrimesh<> *mesh) override {}

  explicit TrettnerSolver(const std::string &path) : mesh_path(path) {}

  void preprocess() override {}

  void query(const int vid, std::vector<double> &res) override {
    res.clear();
    vector<int> vids = {vid};
    half_edge = HEInit(mesh_path, vids);
    res = distance_field_trettner(half_edge, vids, time);
    cout << "? Trettner RES size: " << res.size() << endl;

    // for (int i = 0; i < res.size(); i++) {
    //   cout << "vertex: " << i << ", value: " << res[i] << endl;
    // }
  }
};

// ---------- FAST MARCHING ----------
class FastMarchingSolver : public GeodesicMethod {
public:
  FastMarchingSolver() {}
  ~FastMarchingSolver() {}

  double time;
  vector<double> coords;
  vector<uint> tris;
  string mesh_path;

  // MATLAB engine members
  std::unique_ptr<MATLABEngine> matlabPtr;
  ArrayFactory factory;
  std::unique_ptr<TypedArray<double>> verticesMat;
  std::unique_ptr<TypedArray<double>> facesMat;
  std::unique_ptr<Array> startPointsMat;
  std::unique_ptr<StructArray> options; // Use smart pointer for StructArray
  std::vector<Array> results;

  explicit FastMarchingSolver(const std::string &path) : mesh_path(path) {}

  void load(DrawableTrimesh<> *m) override {
    this->coords = extract_coords(*m);
    this->tris = extract_tris(*m);
  }

  void preprocess() override {

    size_t V = coords.size() / 3;
    size_t F = tris.size() / 3;
    verticesMat = std::make_unique<TypedArray<double>>(
        factory.createArray<double>({V, 3}));
    for (size_t i = 0; i < V; ++i) {
      (*verticesMat)[i][0] = coords[i * 3 + 0];
      (*verticesMat)[i][1] = coords[i * 3 + 1];
      (*verticesMat)[i][2] = coords[i * 3 + 2];
    }

    facesMat = std::make_unique<TypedArray<double>>(
        factory.createArray<double>({F, 3}));
    for (size_t i = 0; i < F; ++i) {
      (*facesMat)[i][0] = static_cast<double>(tris[i * 3 + 0] + 1);
      (*facesMat)[i][1] = static_cast<double>(tris[i * 3 + 1] + 1);
      (*facesMat)[i][2] = static_cast<double>(tris[i * 3 + 2] + 1);
    }

    options = std::make_unique<StructArray>(
        factory.createStructArray({1, 1}, {"verbose"}));
    (*options)[0]["verbose"] = factory.createScalar<double>(1);
  }

  void query(const int vid, std::vector<double> &res) override {
    res.clear();
    time = 0.0;
    matlabPtr = startMATLAB();
    matlabPtr->eval(u"addpath('/Users/giuliabenvenuto/Library/Application "
                    u"Support/MathWorks/MATLAB Add-Ons/Collections/Toolbox "
                    u"Fast Marching/toolbox_fast_marching');");

    startPointsMat =
        std::make_unique<Array>(factory.createScalar<double>(vid + 1));

    auto start = std::chrono::high_resolution_clock::now();
    results =
        matlabPtr->feval(u"perform_fast_marching_mesh", 3,
                         {*verticesMat, *facesMat, *startPointsMat, *options});
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    // save elapsed time
    time = elapsed.count();

    TypedArray<double> D = results[0];
    res = vector<double>(D.begin(), D.end());

    cout << "? FastMarching RES size: " << res.size() << endl;
    cout << "Computation Time FAST MARCHING: " << elapsed.count() << " s"
         << endl;

    // for (int i = 0; i < res.size(); i++) {
    //   cout << "vertex: " << i << ", value: " << res[i] << endl;
    // }
  }

  // function to get the elapsed time
  double get_time() { return time; }
};

// ---------- Heat GEOMETRY CENTRAL ----------
class HeatSolver : public GeodesicMethod {
public:
  HeatSolver() {}
  ~HeatSolver() {}

  gc_mesh gc_m;
  unique_ptr<HeatMethodDistanceSolver> heatSolverGC;
  double time_scalar = 2.0;

  void load(DrawableTrimesh<> *mesh) override {
    gc_m = make_gc_mesh(extract_tris(*mesh), mesh->vector_verts());
  }

  void preprocess() override {
    // Ensure that the mesh and geometry are loaded and valid
    if (!gc_m.topology || !gc_m.geometry) {
      cerr << "Mesh or geometry not initialized." << endl;
      return;
    }

    // Initialize the heat method distance solver
    // heatSolverGC = make_unique<HeatMethodDistanceSolver>(*gc_m.geometry);
    // cout << "Heat time in preprocess: " << time_scalar << endl;
    heatSolverGC =
        make_unique<HeatMethodDistanceSolver>(*gc_m.geometry, time_scalar);
  }

  void set_t(const float new_t) {
    time_scalar = new_t;
    // Create the solver considering the new time scalar
    cout << "Time scalar in set_t: " << time_scalar << endl;
    heatSolverGC =
        make_unique<HeatMethodDistanceSolver>(*gc_m.geometry, time_scalar);
  }

  void query(const int vid, std::vector<double> &res) override {
    cout << "HEAT TIME SCALAR: " << time_scalar << endl;

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

    Vertex sourceVertex =
        Vertex(gc_m.topology.get(), vid); // Convert int vid to Vertex
    VertexData<double> distances = heatSolverGC->computeDistance(sourceVertex);

    res.clear();
    res.resize(distances.size());
    for (size_t i = 0; i < distances.size(); ++i) {
      res[i] = distances[Vertex(gc_m.topology.get(), i)];
    }
    cout << "? Heat RES size: " << res.size() << endl;

    // // print the results
    // for (int i = 0; i < res.size(); i++) {
    //   cout << "vertex: " << i << ", value: " << res[i] << endl;
    // }
  }
};

// ---------- FMM GEOMETRY CENTRAL ----------
class FastMarchingGC : public GeodesicMethod {
public:
  FastMarchingGC() {}
  ~FastMarchingGC() {}

  gc_mesh gc_m;

  void load(DrawableTrimesh<> *mesh) override {
    gc_m = make_gc_mesh(extract_tris(*mesh), mesh->vector_verts());
  }

  void preprocess() override {}

  void query(const int vid, std::vector<double> &res) override {

    // Check if the vertex ID is valid
    if (vid < 0 || vid >= gc_m.topology->nVertices()) {
      cerr << "Invalid vertex ID." << endl;
      return;
    }

    Vertex sourceVertex =
        Vertex(gc_m.topology.get(), vid); // Convert int vid to Vertex
    std::vector<std::pair<Vertex, double>> initialDistances = {
        make_pair(sourceVertex, 0.0)};
    VertexData<double> distances =
        FMMDistance(*gc_m.geometry, initialDistances);

    res.clear();
    res.resize(distances.size());
    for (size_t i = 0; i < distances.size(); ++i) {
      res[i] = distances[Vertex(gc_m.topology.get(), i)];
    }
    // cout << "? FMM RES size: " << res.size() << endl;

    // // print the results
    // for (int i = 0; i < res.size(); i++) {
    //   cout << "vertex: " << i << ", value: " << res[i] << endl;
    // }
  }
};

// ---------- Geotangle ----------
class GeotangleSolver : public GeodesicMethod {
public:
  GeotangleSolver() {}
  ~GeotangleSolver() {}

  DrawableTrimesh<> *m;
  geodesic_solver solver;
  bool solver_computed = false;

  void load(DrawableTrimesh<> *mesh) override {
    m = mesh;
    cout << "Number of vertices: " << m->num_verts() << endl;
    cout << "Number of faces: " << m->num_polys() << endl;
  }

  void preprocess() override {
    solver = make_geodesic_solver(*m, true);
    solver_computed = true;
  }

  void query(const int vid, std::vector<double> &res) override {
    res.clear();
    res = compute_geodesic_distances(solver, {vid});

    // for (int i = 0; i < res.size(); i++) {
    //   cout << "vertex: " << i << ", value: " << res[i] << endl;
    // }
  }
};

// ---------- Edge ----------
class EdgeSolver : public GeodesicMethod {
public:
  EdgeSolver() {}
  ~EdgeSolver() {}

  DrawableTrimesh<> *m;
  geodesic_solver solver;
  bool solver_computed = false;

  void load(DrawableTrimesh<> *mesh) override { m = mesh; }
  void preprocess() override {
    solver = make_geodesic_solver(*m, false);
    solver_computed = true;
  }

  void query(const int vid, std::vector<double> &res) override {
    res.clear();
    res = compute_geodesic_distances(solver, {vid});
    cout << "? Edge RES size: " << res.size() << endl;

    // for (int i = 0; i < res.size(); i++) {
    //   cout << "vertex: " << i << ", value: " << res[i] << endl;
    // }

    // for (int i = 0; i < res.size(); i++) {
    //   cout << i << "," << res[i] << endl;
    // }
  }
};

// ---------- Extended ----------
class ExtendedSolver : public GeodesicMethod {
public:
  ExtendedSolver() {}
  ~ExtendedSolver() {}

  DrawableTrimesh<> *m;
  geodesic_solver solver;
  dual_geodesic_solver dual_solver;
  bool dual_solver_computed = false;
  int k = 3;

  void load(DrawableTrimesh<> *mesh) override { m = mesh; }

  void set_k(const int new_k, const bool compute_solver = true) {
    k = new_k;
    if (compute_solver) {
      if (!dual_solver_computed) {
        dual_solver = make_dual_geodesic_solver(*m);
        dual_solver_computed = true;
      }
      solver = extended_solver(*m, dual_solver, k);
    }
  }

  void preprocess() override {
    cout << "K: " << k << endl;
    cout << "Make dual geodesic solver START" << endl;
    dual_solver = make_dual_geodesic_solver(*m);
    cout << "Make dual geodesic solver END" << endl;

    cout << "Extended solver START" << endl;
    solver = extended_solver(*m, dual_solver, k);
    cout << "Extended solver END" << endl;
    dual_solver_computed = true;
  }

  void query(const int vid, std::vector<double> &res) override {
    res.clear();
    res = compute_geodesic_distances(solver, {vid});

    // for (int i = 0; i < res.size(); i++) {
    //   cout<< i << "," << res[i] << endl;
    // }

    cout << "? Geotangle RES size: " << res.size() << endl;
  }
};

// ---------- Lanthier ----------
class LanthierSolver : public GeodesicMethod {
public:
  LanthierSolver() {}
  ~LanthierSolver() {}

  DrawableTrimesh<> *m;
  geodesic_solver solver;
  int n_steiner = 9;

  void load(DrawableTrimesh<> *mesh) override { m = mesh; }

  void set_n_steiner(const int new_n_steiner) {
    n_steiner = new_n_steiner;
    solver = compute_fine_graph(*m, n_steiner);
  }

  void preprocess() override { solver = compute_fine_graph(*m, n_steiner); }

  void query(const int vid, std::vector<double> &res) override {
    res.clear();
    res = compute_geodesic_distances(solver, {vid});
    // for (int i = 0; i < res.size(); i++) {
    //   cout << "vertex: " << i << ", value: " << res[i] << endl;
    // }
    res.resize(m->num_verts());
    cout << "? Lanthier RES size: " << res.size() << endl;
  }
};

#endif