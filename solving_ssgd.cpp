// Old implementation of the SSGD methods
#include "solving_ssgd.h"
#include "SSGD_methods/Graph-based_methods/shortest_path.h"


ScalarField SSGD_Heat(DrawableTrimesh<> &m,
                      GeodesicsCache &prefactored_matrices,
                      vector<uint> &sources, double &time_heat) {

  bool cache = false;
  if (prefactored_matrices.heat_flow_cache != NULL) {
    cache = true;
  }
  ScalarField sc_heat;
  // Timer
  auto start_heat = chrono::high_resolution_clock::now();
  // Method inside: #include <cinolib/geodesics.h>
  sc_heat = compute_geodesics_amortized(m, prefactored_matrices, sources);
  auto stop_heat = chrono::high_resolution_clock::now();
  time_heat =
      chrono::duration_cast<chrono::milliseconds>(stop_heat - start_heat)
          .count();

  if (cache) {
    cout << "Heat computation with cache: " << time_heat << " milliseconds"
         << endl;
  } else {
    cout << "Heat computation without cache: " << time_heat << " milliseconds"
         << endl;
  }

  // return compute_geodesics_amortized(m, prefactored_matrices, sources);
  return sc_heat;
}

// SSGD_VTP method with multiple sources
ScalarField SSGD_VTP(DrawableTrimesh<> &m, vector<int> &sources,
                     double &vtp_geodesic_time) {

  vector<double> field_data(m.num_verts(), std::numeric_limits<double>::max());
  ScalarField sc_vtp;

  // Only one source
  // auto start_geodesic_VTP = chrono::high_resolution_clock::now();
  // field_data = exact_geodesic_distance(m.vector_polys(), m.vector_verts(),
  // sources[0]); auto stop_geodesic_VTP = chrono::high_resolution_clock::now();

  auto start_geodesic_VTP = chrono::high_resolution_clock::now();
  // Iterate over each source
  for (int source : sources) {
    vector<double> current_distances = exact_geodesic_distance(m.vector_polys(), m.vector_verts(), source);

    // Update the minimum distance for each vertex
    for (size_t i = 0; i < field_data.size(); ++i) {
      field_data[i] = std::min(field_data[i], current_distances[i]);
    }
  }

  auto stop_geodesic_VTP = chrono::high_resolution_clock::now();

  // Invert the color mapping
  for (auto &value : field_data) {
    value = 1.0 - value;
  }

  sc_vtp = ScalarField(field_data);
  sc_vtp.normalize_in_01();

  vtp_geodesic_time = chrono::duration_cast<chrono::milliseconds>(
                          stop_geodesic_VTP - start_geodesic_VTP)
                          .count();
  cout << "Geodesic computation with VTP: " << vtp_geodesic_time
       << " milliseconds" << endl;

  return sc_vtp;
}

ScalarField SSGD_GeoTangle(DrawableTrimesh<> &m, geodesic_solver &solver,
                           vector<int> &sources,
                           double &geotangle_geodesic_time) {

  vector<double> distances_geo;
  ScalarField sc_geotangle;

  auto start_geodesic_GeoTangle = chrono::high_resolution_clock::now();
  distances_geo = compute_geodesic_distances(solver, sources);
  auto stop_geodesic_GeoTangle = chrono::high_resolution_clock::now();

  // Invert the color mapping
  for (auto &value : distances_geo) {
    value = 1.0 - value;
  }

  sc_geotangle = ScalarField(distances_geo);
  sc_geotangle.normalize_in_01();

  auto duration_geodesic_GeoTangle =
      chrono::duration_cast<chrono::milliseconds>(stop_geodesic_GeoTangle -
                                                  start_geodesic_GeoTangle);
  geotangle_geodesic_time =
      chrono::duration_cast<chrono::milliseconds>(stop_geodesic_GeoTangle -
                                                  start_geodesic_GeoTangle)
          .count();
  cout << "Geodesic computation with GeoTangle: "
       << duration_geodesic_GeoTangle.count() << " milliseconds" << endl;

  return sc_geotangle;
}

ScalarField SSGD_Edge(DrawableTrimesh<> &m, geodesic_solver &solver,
                      vector<int> &sources, double &edge_geodesic_time) {

  vector<double> distances_edge;
  ScalarField sc_edge;

  auto start_geodesic_edge = chrono::high_resolution_clock::now();
  distances_edge = compute_geodesic_distances(solver, sources);
  auto stop_geodesic_edge = chrono::high_resolution_clock::now();

  // Invert the color mapping
  for (auto &value : distances_edge) {
    value = 1.0 - value;
  }

  sc_edge = ScalarField(distances_edge);
  sc_edge.normalize_in_01();

  auto duration_geodesic_edge = chrono::duration_cast<chrono::milliseconds>(
      stop_geodesic_edge - start_geodesic_edge);
  edge_geodesic_time = chrono::duration_cast<chrono::milliseconds>(
                           stop_geodesic_edge - start_geodesic_edge)
                           .count();
  cout << "Geodesic computation with Edge: " << duration_geodesic_edge.count()
       << " milliseconds" << endl;

  return sc_edge;
}

ScalarField SSGD_Extended(DrawableTrimesh<> &m, geodesic_solver &solver,
                          vector<int> &sources,
                          double &extended_geodesic_time) {

  vector<double> distances_ext(
      m.num_verts(),
      std::numeric_limits<double>::max()); // Initialize with max double value
  ScalarField sc_extended;

  auto start_geodesic_extended = chrono::high_resolution_clock::now();
  distances_ext = compute_geodesic_distances(solver, sources);
  auto stop_geodesic_extended = chrono::high_resolution_clock::now();

  // Invert the color mapping
  for (auto &value : distances_ext) {
    value = 1.0 - value;
  }

  sc_extended = ScalarField(distances_ext);
  sc_extended.normalize_in_01();

  auto duration_geodesic_extended = chrono::duration_cast<chrono::milliseconds>(
      stop_geodesic_extended - start_geodesic_extended);
  extended_geodesic_time = chrono::duration_cast<chrono::milliseconds>(
                               stop_geodesic_extended - start_geodesic_extended)
                               .count();
  cout << "Geodesic computation with Extended: "
       << duration_geodesic_extended.count() << " milliseconds" << endl;

  return sc_extended;
}
