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
    cout << "VTP query started... " << endl;
    res.clear();  // Clear previous results

    auto start = std::chrono::high_resolution_clock::now();
    res = exact_geodesic_distance(m.vector_polys(), m.vector_verts(), vid);
    auto end = std::chrono::high_resolution_clock::now();
    // print time in seconds with 4 decimal digits
    std::chrono::duration<double> elapsed = end - start;
    cout << "VTP computation: " << elapsed.count() << " s" << endl;



    auto start_modify = std::chrono::high_resolution_clock::now();
    for (auto &value : res) {
      value = 1.0 - value;
    }
    auto end_modify = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_modify = end_modify - start_modify;
    cout << "Modify the results: " << elapsed_modify.count() << " s" << endl;

    auto start_scalar = std::chrono::high_resolution_clock::now();
    sc = ScalarField(res);
    sc.normalize_in_01();
    auto end_scalar = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_scalar = end_scalar - start_scalar;
    cout << "ScalarField computation: " << elapsed_scalar.count() << " s" << endl;
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

    res.clear();  // Clear previous results
    res.resize(sc.size());  // Ensure 'res' can hold all distances
    for (size_t i = 0; i < sc.size(); ++i) {
        res[i] = sc[i];  // Copy distances to 'res'
    }
  }

};


// ---------- FAST MARCHING ----------
class FastMarchingSolver : public GeodesicMethod {
public:
  FastMarchingSolver() {}
  ~FastMarchingSolver() {}

  DrawableTrimesh<> m;
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

  void load(const std::vector<double> &coords, const std::vector<uint> &tris) override {
    m = DrawableTrimesh(coords, tris);
    this->coords = coords;
    this->tris = tris;
  }

  void preprocess() override {

    verticesMat = std::make_unique<TypedArray<double>>(factory.createArray<double>({m.num_verts(), 3}));
    for (size_t i = 0; i < m.num_verts(); ++i) {
        (*verticesMat)[i][0] = coords[i * 3 + 0];
        (*verticesMat)[i][1] = coords[i * 3 + 1];
        (*verticesMat)[i][2] = coords[i * 3 + 2];
    }

    facesMat = std::make_unique<TypedArray<double>>(factory.createArray<double>({m.num_polys(), 3}));
    for (size_t i = 0; i < m.num_polys(); ++i) {
        (*facesMat)[i][0] = static_cast<double>(tris[i * 3 + 0] + 1);
        (*facesMat)[i][1] = static_cast<double>(tris[i * 3 + 1] + 1);
        (*facesMat)[i][2] = static_cast<double>(tris[i * 3 + 2] + 1);
    }

    options = std::make_unique<StructArray>(factory.createStructArray({1, 1}, {"verbose"}));
    (*options)[0]["verbose"] = factory.createScalar<double>(1);
  }

  void query(const int vid, std::vector<double> &res, ScalarField &sc) override {
    matlabPtr = startMATLAB();
    matlabPtr->eval(u"addpath('/Users/giuliabenvenuto/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/Toolbox Fast Marching/toolbox_fast_marching');");

    
    startPointsMat = std::make_unique<Array>(factory.createScalar<double>(vid + 1));
    results = matlabPtr->feval(u"perform_fast_marching_mesh", 3, {*verticesMat, *facesMat, *startPointsMat, *options});
    TypedArray<double> D = results[0];
    std::vector<double> distanceField(D.begin(), D.end());
    for (auto &value : distanceField) {
        value = 1.0 - value; // Invert the value
    }
    sc = ScalarField(distanceField);
    sc.normalize_in_01();

    res.clear();  // Clear previous results
    res.resize(sc.size());  // Ensure 'res' can hold all distances
    for (size_t i = 0; i < sc.size(); ++i) {
        res[i] = sc[i];  // Copy distances to 'res'
    }
  }
};


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
    res.clear(); 
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
    res.clear();
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
    // cout << "K in the function: " << k << endl;
    if (compute_solver) {
      if (!dual_solver_computed) {
        dual_solver = make_dual_geodesic_solver(m);
        dual_solver_computed = true;
      }
      solver = extended_solver(m, dual_solver, k);
    }
  }

  void preprocess() override {
    cout << "Make dual geodesic solver START" << endl;
    dual_solver = make_dual_geodesic_solver(m); // TODO: risolvi che questo blocca delle mesh
    cout << "Make dual geodesic solver END" << endl;

    cout << "Extended solver START" << endl;
    solver = extended_solver(m, dual_solver, k);
    cout << "Extended solver END" << endl;
    dual_solver_computed = true;
  }

  void query(const int vid, std::vector<double> &res, ScalarField &sc) override {
    res.clear();
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
  geodesic_solver solver;
  //dual_geodesic_solver dual_solver;
  //bool dual_solver_computed = false;
  int n_steiner = 1;

  void load(const std::vector<double> &coords, const std::vector<uint> &tris) override {
    m = DrawableTrimesh(coords, tris);
  }
  
  void set_n_steiner(const int new_n_steiner) {
    n_steiner = new_n_steiner;
    //cout << "n_steiner in the function: " << n_steiner << endl;
    solver = compute_fine_graph(m, n_steiner);

  }

  void preprocess() override {
    solver = compute_fine_graph(m, n_steiner);
  }

  void query(const int vid, std::vector<double> &res, ScalarField &sc) override {
    // solver computed with the lanthier method
    // res -> vector<double>
    res.clear();
    res = compute_geodesic_distances(solver, {vid});
    // get only m.num_verts() elements of res
    res.resize(m.num_verts());
    for (auto &value : res) {
      value = 1.0 - value;
    }
    sc = ScalarField(res);
    sc.normalize_in_01();
  }

};

#endif