#ifndef GEO_METHOD
#define GEO_METHOD
#include <cinolib/how_many_seconds.h>
#include <functional>
#include <sys/types.h>

using namespace cinolib;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

////// ABSTRACT CLASS //////
class GeodesicMethod {
public:
  explicit GeodesicMethod() {}
  virtual ~GeodesicMethod() {}

  virtual void load(const std::vector<double> &coords,
                    const std::vector<uint> &tris) = 0;
  virtual void preprocess() = 0;
  virtual void query(const uint vid, std::vector<double> &res) = 0;
};

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

////// INSTANCE OF A SPECIFIC ALGORITHM //////
class MyMethod : public GeodesicMethod {
public:
  MyMethod() {}
  ~MyMethod() {}

  void load(const std::vector<double> &coords,
            const std::vector<uint> &tris) override {
    // mesh loading
  }

  void preprocess() override {
    // preprocessing
  }

  void query(const uint vid, std::vector<double> &res) override {
    // query
  }
};

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

////// BENCHMARK FUNCTION THAT OPERATES ON A GENERIC (UNKNOWN) ALGORITHM //////
void benchmark(GeodesicMethod &m) {
  std::vector<double> mesh_coords;
  std::vector<uint> mesh_tris;
  std::chrono::steady_clock::time_point tic, toc;

  tic = std::chrono::steady_clock::now();
  m.load(mesh_coords, mesh_tris);
  toc = std::chrono::steady_clock::now();
  double t_load = how_many_seconds(tic, toc);

  tic = std::chrono::steady_clock::now();
  m.preprocess();
  toc = std::chrono::steady_clock::now();
  double t_pre = how_many_seconds(tic, toc);

  std::vector<double> res;
  tic = std::chrono::steady_clock::now();
  m.query(0, res);
  toc = std::chrono::steady_clock::now();
  double t_query = how_many_seconds(tic, toc);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#endif