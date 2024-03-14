#include "solving_ssgd.h"


ScalarField SSGD_Heat(DrawableTrimesh<> &m, GeodesicsCache &prefactored_matrices, vector<uint> &sources, double &time_heat) {
  
    bool cache = false;
    if (prefactored_matrices.heat_flow_cache != NULL) {
    cache = true;
    }
    ScalarField sc_heat;
    // Timer
    auto start_heat = chrono::high_resolution_clock::now();
    // Method inside: #include <cinolib/geodesics.h>
    // compute_geodesics_amortized(m, prefactored_matrices, sources).copy_to_mesh(m);
    sc_heat = compute_geodesics_amortized(m, prefactored_matrices, sources);
    auto stop_heat = chrono::high_resolution_clock::now();

    time_heat = chrono::duration_cast<chrono::milliseconds>(stop_heat - start_heat).count();

    // m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);

    if (cache) {
    cout << "Heat computation with cache: " << time_heat << " milliseconds" << endl;
    } else {
    cout << "Heat computation without cache: " << time_heat << " milliseconds" << endl;
    }

    //return compute_geodesics_amortized(m, prefactored_matrices, sources);
    return sc_heat;
}


ScalarField SSGD_VTP(DrawableTrimesh<> &m, vector<int> &sources, double &vtp_geodesic_time) {
    vector<double> field_data;
    ScalarField sc_vtp;

    auto start_geodesic_VTP = chrono::high_resolution_clock::now();
    field_data = exact_geodesic_distance(m.vector_polys(), m.vector_verts(), sources[0]);
    auto stop_geodesic_VTP = chrono::high_resolution_clock::now();

    // Invert the color mapping
    for (auto &value : field_data) {
        value = 1.0 - value;
    }

    sc_vtp = ScalarField(field_data);
    sc_vtp.normalize_in_01();
    //field.copy_to_mesh(m);
    //m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);

    auto duration_geodesic_VTP = chrono::duration_cast<chrono::milliseconds>(stop_geodesic_VTP - start_geodesic_VTP);
    vtp_geodesic_time = chrono::duration_cast<chrono::milliseconds>(stop_geodesic_VTP - start_geodesic_VTP).count();
    cout << "Geodesic computation with VTP: " << duration_geodesic_VTP.count()
        << " milliseconds" << endl;

    return sc_vtp;
} 


ScalarField SSGD_GeoTangle(DrawableTrimesh<> &m, geodesic_solver &solver, vector<int> &sources, double &geotangle_geodesic_time) {

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
  //field_geo.copy_to_mesh(m);
  //m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);

  auto duration_geodesic_GeoTangle = chrono::duration_cast<chrono::milliseconds>(stop_geodesic_GeoTangle - start_geodesic_GeoTangle);
  geotangle_geodesic_time = chrono::duration_cast<chrono::milliseconds>(stop_geodesic_GeoTangle - start_geodesic_GeoTangle).count();
  cout << "Geodesic computation with GeoTangle: " << duration_geodesic_GeoTangle.count() << " milliseconds" << endl;

  return sc_geotangle;
}


ScalarField SSGD_Edge(DrawableTrimesh<> &m, geodesic_solver &solver, vector<int> &sources, double &edge_geodesic_time) {

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
  //field_edge.copy_to_mesh(m);
  //m.show_texture1D(TEXTURE_1D_HSV_W_ISOLINES);

  auto duration_geodesic_edge = chrono::duration_cast<chrono::milliseconds>(stop_geodesic_edge - start_geodesic_edge);
  edge_geodesic_time = chrono::duration_cast<chrono::milliseconds>(stop_geodesic_edge - start_geodesic_edge).count();
  cout << "Geodesic computation with Edge: " << duration_geodesic_edge.count()
       << " milliseconds" << endl;

    return sc_edge;
}