#include "diff_geo.h"

using namespace Eigen;

vec3d project_vec(const vec3d &v, const vec3d &n) {
  auto proj = n * v.dot(n);

  return v - proj;
}
vec3d polar_x_axis(const DrawableTrimesh<> &m, const int vid, const vec3d &N) {
  uint vid0 = m.adj_v2v(vid)[0];
  vec3d v = m.vert(vid0) - m.vert(vid);
  vec3d e = project_vec(v, N);
  e.normalize();
  return e;
}
vec3d average_normal(const vector<vec3d> &normals, const vector<int> &nbr,
                     const int vid) {
  vec3d result = vec3d(0, 0, 0);
  for (int nei : nbr) {
    result += normals[nei];
  }
  return result / (double)nbr.size();
}
vector<int> k_ring(const DrawableTrimesh<> &m, const int vid, const int k) {

  vector<int> active_set = {vid};
  vector<int> ring = active_set;
  for (uint i = 0; i < k; ++i) {
    vector<int> next_active_set;
    for (uint curr : active_set)
      for (uint nbr : m.adj_v2v(curr)) {
        if (find(ring.begin(), ring.end(), nbr) == ring.end() && nbr != vid) {
          next_active_set.push_back(nbr);
          ring.push_back(nbr);
        }
      }

    active_set = next_active_set;
  }

  return ring;
}

std::tuple<vector<int>, vector<double>, vector<double>, vec3d>
filtered_ring_stencil(const DrawableTrimesh<> &m, const int vid,
                      const int k = 1) {

  vector<int> nbr = k_ring(m, vid, k);

  vec3d N = average_normal(m.vector_vert_normals(), nbr, vid);
  vec3d pos = m.vert(vid);
  auto count = 1;
  while (nbr.size() < 5) {
    nbr = k_ring(m, vid, k + count);
    ++count;
  }
  vector<double> lens(nbr.size());
  vector<double> thetas(nbr.size());

  auto e = polar_x_axis(m, vid, N);
  for (size_t i = 0; i < nbr.size(); ++i) {
    vec3d v = m.vert(nbr[i]) - pos;
    auto proj = project_vec(v, N);
    lens[i] = proj.norm();
    double theta = e.angle_rad(proj);
    if (N.dot(e.cross(proj)) < 0)
      theta = 2 * M_PI - theta;
    thetas[i] = (nbr[i] == vid) ? 0 : theta;
  }

  return {nbr, lens, thetas, N};
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// Patch fitting
vec3d evaluate_quadric(const MatrixXd &Q, const vec2d &uv) {

  auto pos = vec3d{0, 0, 0};
  for (auto j = 0; j < 3; ++j) {
    pos[j] = Q(0, j) + uv.x() * Q(1, j) + uv.y() * Q(2, j) +
             1.f / 2 * Q(3, j) * std::pow(uv.x(), 2) +
             Q(4, j) * uv.x() * uv.y() +
             1.f / 2 * Q(5, j) * std::pow(uv.y(), 2);
  }

  return pos;
}

vec3d evaluate_quadric_du(const MatrixXd &Q, const vec2d &uv) {

  auto pos = vec3d{0, 0, 0};
  for (auto j = 0; j < 3; ++j) {
    pos[j] = Q(1, j) + Q(3, j) * uv.x() + Q(4, j) * uv.y();
  }

  return pos;
}
vec3d evaluate_quadric_dv(const MatrixXd &Q, const vec2d &uv) {

  auto pos = vec3d{0, 0, 0};
  for (auto j = 0; j < 3; ++j) {
    pos[j] = Q(2, j) + Q(5, j) * uv.y() + Q(4, j) * uv.x();
  }

  return pos;
}
vec3d evaluate_quadric_duu(const MatrixXd &Q, const vec2d &uv) {

  auto pos = vec3d{0, 0, 0};
  for (auto j = 0; j < 3; ++j) {
    pos[j] = Q(3, j);
  }

  return pos;
}
vec3d evaluate_quadric_dvv(const MatrixXd &Q, const vec2d &uv) {

  auto pos = vec3d{0, 0, 0};
  for (auto j = 0; j < 3; ++j) {
    pos[j] = Q(5, j);
  }

  return pos;
}
vec3d evaluate_quadric_duv(const MatrixXd &Q, const vec2d &uv) {

  auto pos = vec3d{0, 0, 0};
  for (auto j = 0; j < 3; ++j) {
    pos[j] = Q(4, j);
  }

  return pos;
}
Matrix2d first_fund(const patch &p) {
  Matrix2d I;

  I << p.xu.dot(p.xu), p.xu.dot(p.xv), p.xu.dot(p.xv), p.xv.dot(p.xv);

  return I;
}
Matrix2d second_fund(const patch &p) {
  Matrix2d II;

  vec3d N = p.xu.cross(p.xv);
  N.normalize();
  II << N.dot(p.xuu), N.dot(p.xuv), N.dot(p.xuv), N.dot(p.xvv);

  return II;
}
Matrix2d shape_operator(const patch &p) {
  auto I = first_fund(p);
  auto II = second_fund(p);
  assert(I.determinant() > 1e-10);
  Matrix2d A = -I.inverse() * II;

  return A;
}
std::tuple<vec2d, vector<vec2d>> principal_curv_and_dir(const patch &p) {
  auto S = shape_operator(p);
  SelfAdjointEigenSolver<Matrix2d> dec(S);
  auto eivals = dec.eigenvalues();
  auto eivec = dec.eigenvectors();
  double k1 = eivals(0);
  double k2 = eivals(1);

  Vector2d K1 = eivec.col(0);
  Vector2d K2 = eivec.col(1);
  auto d1 = vec2d{K1(0), K1(1)};
  auto d2 = vec2d{K2(0), K2(1)};
  d1.normalize();
  d2.normalize();
  return {vec2d{k1, k2}, {d1, d2}};
}
std::tuple<vector<vec3d>, vector<vec3d>>
principal_curvature_field(const vector<patch> &quadrics) {
  auto n = quadrics.size();
  auto k1 = vector<vec3d>(n);
  auto k2 = vector<vec3d>(n);
  for (auto i = 0; i < n; ++i) {
    auto curr_patch = quadrics[i];
    auto xu = evaluate_quadric_du(curr_patch.quadric, vec2d(0, 0));
    auto xv = evaluate_quadric_dv(curr_patch.quadric, vec2d(0, 0));
    auto [k, d] = principal_curv_and_dir(curr_patch);
    k1[i] = d[0].x() * xu + d[0].y() * xv;
    k2[i] = d[1].x() * xu + d[1].y() * xv;
  }
  return {k1, k2};
}

Matrix2d isophotic_metric(const patch &p, const double &w = 0,
                          const double &w_star = 1) {
  Matrix2d I;
  Matrix2d II;
  vec3d n = p.xu.cross(p.xv);
  n.normalize();
  double E = p.xu.dot(p.xu);
  double F = p.xu.dot(p.xv);
  double G = p.xv.dot(p.xv);
  double L = n.dot(p.xuu);
  double M = n.dot(p.xuv);
  double N = n.dot(p.xvv);
  I << E, F, F, G;
  II << L, M, M, N;

  double K = (L * N - pow(M, 2)) / (E * G - pow(F, 2));
  double H =
      (E * N - 2 * F * M + G * L) / (E * G - pow(F, 2)); // multplied by 2

  Matrix2d III = H * II - K * I;

  return w * I + w_star * III;
}
double max_residual(const MatrixXd &r) {
  auto max_err = __DBL_MIN__;
  for (auto i = 0; i < r.rows(); ++i) {
    vec3d curr_err = vec3d(r(i, 0), r(i, 1), r(i, 2));
    max_err = std::max(curr_err.norm(), max_err);
  }
  return max_err;
}
patch patch_fitting(const DrawableTrimesh<> &m, const int vid, const int k) {

  patch result;
  double u_max, v_max, u_min, v_min, w_max, w_min;
  u_max = v_max = w_max = __DBL_MIN__;
  u_min = v_min = w_min = __DBL_MAX__;
  auto [nbr, lens, tetas, n] = filtered_ring_stencil(m, vid, k);
  size_t s = nbr.size();
  vec3d vert = m.vert(vid);
  Eigen::MatrixXd A(s, 6);
  Eigen::MatrixXd P(s, 3);
  auto &p_nbr = result.parametric_nbr;
  p_nbr.resize(s);
  for (size_t i = 0; i < s; ++i) {
    vec2d pos = vec2d{lens[i] * cos(tetas[i]), lens[i] * sin(tetas[i])};
    u_max = std::max(u_max, pos.x());
    v_max = std::max(v_max, pos.y());
    u_min = std::min(u_min, pos.x());
    v_min = std::min(v_min, pos.y());
    vec3d coords = m.vert(nbr[i]);
    p_nbr[i] = make_pair(pos, nbr[i]);
    for (auto h = 0; h < 3; ++h)
      P(i, h) = coords[h];

    A(i, 0) = 1;
    A(i, 1) = pos[0];
    A(i, 2) = pos[1];
    A(i, 3) = 1.f / 2 * pow(pos[0], 2);
    A(i, 4) = pos[0] * pos[1];
    A(i, 5) = 1.f / 2 * pow(pos[1], 2);
  }

  MatrixXd At = Transpose<MatrixXd>(A);
  MatrixXd B = At * A;
  ColPivHouseholderQR<MatrixXd> dec(B);
  MatrixXd c(6, 3);

  for (auto h = 0; h < 3; ++h) {
    c.col(h) = dec.solve(At * P.col(h));
  }

  MatrixXd res = A * c - P;
  result.max_res = max_residual(res);
  result.xu = evaluate_quadric_du(c, vec2d(0, 0));
  result.xv = evaluate_quadric_dv(c, vec2d(0, 0));
  result.xuu = evaluate_quadric_duu(c, vec2d(0, 0));
  result.xuv = evaluate_quadric_duv(c, vec2d(0, 0));
  result.xvv = evaluate_quadric_dvv(c, vec2d(0, 0));
  result.domain_range = {vec2d{u_min, u_max}, vec2d{v_min, v_max}};
  result.quadric = c;

  return result;
}
vector<patch> patch_fitting(const DrawableTrimesh<> &m, const int k) {

  vector<patch> result(m.num_verts());
  for (uint i = 0; i < m.num_verts(); ++i)
    result[i] = patch_fitting(m, i, k);

  return result;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// Geodesic Solver
template <typename Update, typename Stop, typename Exit>
void visit_geodesic_graph(vector<double> &field, const geodesic_solver &solver,
                          const vector<int> &sources, const int type,
                          Update &&update, Stop &&stop, Exit &&exit) {
  /*
     This algortithm uses the heuristic Small Label Fisrt and Large Label Last
     https://en.wikipedia.org/wiki/Shortest_Path_Faster_Algorithm

     Large Label Last (LLL): When extracting nodes from the queue, pick the
     front one. If it weights more than the average weight of the queue, put
     on the back and check the next node. Continue this way.
     Sometimes average_weight is less than every value due to floating point
     errors (doesn't happen with double precision).

     Small Label First (SLF): When adding a new node to queue, instead of
     always pushing it to the end of the queue, if it weights less than the
     front node of the queue, it is put on front. Otherwise the node is put at
     the end of the queue.
  */

  auto in_queue = vector<bool>(solver.graph.size(), false);

  // Cumulative weights of elements in queue. Used to keep track of the
  // average weight of the queue.
  auto cumulative_weight = 0.0;

  // setup queue
  auto queue = deque<int>();
  for (auto source : sources) {
    in_queue[source] = true;
    cumulative_weight += field[source];
    queue.push_back(source);
  }

  while (!queue.empty()) {
    auto node = queue.front();
    auto average_weight = (float)(cumulative_weight / queue.size());

    // Large Label Last (see comment at the beginning)
    for (auto tries = 0; tries < queue.size() + 1; tries++) {
      if (field[node] <= average_weight)
        break;
      queue.pop_front();
      queue.push_back(node);
      node = queue.front();
    }

    // Remove node from queue.
    queue.pop_front();
    in_queue[node] = false;
    cumulative_weight -= field[node];

    // Check early exit condition.
    if (exit(node))
      break;
    if (stop(node))
      continue;

    for (auto i = 0; i < (int)solver.graph[node].size(); i++) {
      // Distance of neighbor through this node
      double new_distance;
      if (type == geodesic)
        new_distance = field[node] + solver.graph[node][i].length;
      else if (type == isophotic)
        new_distance = field[node] + solver.graph[node][i].isophotic_length;

      auto neighbor = solver.graph[node][i].node;

      auto old_distance = field[neighbor];
      if (new_distance >= old_distance)
        continue;

      if (in_queue[neighbor]) {
        // If neighbor already in queue, don't add it.
        // Just update cumulative weight.
        cumulative_weight += new_distance - old_distance;
      } else {
        // If neighbor not in queue, add node to queue using Small Label
        // First (see comment at the beginning).
        if (queue.empty() || (new_distance < field[queue.front()]))
          queue.push_front(neighbor);
        else
          queue.push_back(neighbor);

        // Update queue information.
        in_queue[neighbor] = true;
        cumulative_weight += new_distance;
      }

      // Update distance of neighbor.
      field[neighbor] = new_distance;
      update(neighbor);
    }
  }
}
void update_geodesic_distances(vector<double> &distances,
                               const geodesic_solver &solver,
                               const vector<int> &sources, const int type,
                               double max_distance = __DBL_MAX__) {

  auto update = [](int node) {};
  auto stop = [&](int node) { return distances[node] > max_distance; };
  auto exit = [](int node) { return false; };
  visit_geodesic_graph(distances, solver, sources, type, update, stop, exit);
}
vector<int> update_geodesic_distances_tracking_propoagation(
    vector<double> &distances, const geodesic_solver &solver,
    const vector<int> &sources, const int type,
    double max_distance = __DBL_MAX__) {

  auto voronoi = vector<int>{};
  auto update = [&voronoi](int node) { voronoi.push_back(node); };
  auto stop = [&](int node) { return distances[node] > max_distance; };
  auto exit = [](int node) { return false; };
  visit_geodesic_graph(distances, solver, sources, type, update, stop, exit);
  return voronoi;
}
vector<double> compute_geodesic_distances(const geodesic_solver &solver,
                                          const vector<int> &sources,
                                          const int type) {

  auto field = vector<double>(solver.graph.size(), __DBL_MAX__);
  for (auto source : sources)
    field[source] = 0.0;

  update_geodesic_distances(field, solver, sources, type);

  return field;
}
std::pair<vector<double>, vector<unsigned>>
exact_geodesic_wrapper(const vector<vector<uint>> &triangles,
                       const vector<vec3d> &positions) {
  int V = (int)positions.size();
  int F = (int)triangles.size();
  vector<double> points(3 * V);
  vector<uint> faces(3 * F);
  vector<double> f(V);

  for (int i = 0; i < V; ++i) {
    points[3 * i] = positions[i].x();
    points[3 * i + 1] = positions[i].y();
    points[3 * i + 2] = positions[i].z();
  }
  for (int i = 0; i < F; ++i) {
    faces[3 * i] = triangles[i][0];
    faces[3 * i + 1] = triangles[i][1];
    faces[3 * i + 2] = triangles[i][2];
  }

  return std::make_pair(points, faces);
}
vector<double> exact_geodesic_distance(const vector<vector<uint>> &triangles,
                                       const vector<vec3d> &positions,
                                       const int &source) {
  int V = (int)positions.size();
  int F = (int)triangles.size();
  vector<double> f;
  f.resize(V);
  auto [points, faces] = exact_geodesic_wrapper(triangles, positions);
  geodesic_VTP::Mesh mesh;
  mesh.initialize_mesh_data(points, faces);
  geodesic_VTP::GeodesicAlgorithmExact algorithm(&mesh);
  algorithm.propagate(source);
  vector<geodesic_VTP::Vertex> verts = mesh.vertices();
  for (int j = 0; j < V; ++j) {
    geodesic_VTP::Vertex v = verts[j];
    auto value = v.geodesic_distance();
    f[j] = value;
  }

  return f;
}
geodesic_solver compute_geodesic_solver(const DrawableTrimesh<> &m,
                                        const vector<patch> &quadrics) {
  geodesic_solver solver;
  uint V = m.num_verts();
  solver.graph.resize(V);
  for (auto i = 0; i < V; ++i) {

    // vector<double> f =
    //     exact_geodesic_distance(m.vector_polys(), m.vector_verts(), i);
    patch p = quadrics[i];
    vec3d pos = m.vert(i);
    auto nbr = p.parametric_nbr;
    size_t k = nbr.size();
    solver.graph[i].resize(k);
    Matrix2d I = first_fund(p);
    Matrix2d S = isophotic_metric(p, 0, 1);
    for (auto j = 0; j < k; ++j) {
      solver.graph[i][j].node = nbr[j].second;
      vec2d v = nbr[j].first;
      Vector2d vd;
      vd << v.x(), v.y();
      solver.graph[i][j].length = vd.transpose() * I * vd;
      solver.graph[i][j].isophotic_length = std::abs(vd.transpose() * S * vd);
    }
  }

  return solver;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// Intrinsic Voronoi diagram
vector<int> update_voronoi_verts(iVd &voronoi, const DrawableTrimesh<> &mesh,
                                 const vector<int> &curr_region) {
  auto updated = vector<int>{};
  auto parsed = vector<bool>(mesh.num_polys(), false);
  for (auto vert : curr_region) {
    auto star = mesh.adj_v2p(vert);
    for (auto tri : star) {
      if (parsed[tri])
        continue;
      parsed[tri] = true;
      auto it = find(voronoi.voronoi_verts.begin(), voronoi.voronoi_verts.end(),
                     (int)tri);

      auto tag_x = voronoi.voronoi_tags[mesh.poly_vert_id(tri, 0)];
      auto tag_y = voronoi.voronoi_tags[mesh.poly_vert_id(tri, 1)];
      auto tag_z = voronoi.voronoi_tags[mesh.poly_vert_id(tri, 2)];

      if (tag_x != tag_y && tag_x != tag_z && tag_y != tag_z) {
        if (it == voronoi.voronoi_verts.end())
          voronoi.voronoi_verts.push_back(tri);
        updated.push_back(tri);
      } else if (it != voronoi.voronoi_verts.end()) {
        voronoi.voronoi_verts.erase(it);
      }
    }
  }
  return updated;
}
void update_voronoi_regions(iVd &vor, const vector<int> &voronoi_centers,
                            const DrawableTrimesh<> &m) {
  vor.voronoi_regions.clear();
  vor.voronoi_regions.resize(voronoi_centers.size());
  for (auto &region : vor.voronoi_regions) {
    region.clear();
    region.reserve(m.num_verts());
  }

  for (auto i = 0; i < m.num_verts(); ++i) {
    auto tag = vor.voronoi_tags[i];
    auto it = find(voronoi_centers.begin(), voronoi_centers.end(), tag);
    if (it == voronoi_centers.end())
      std::cout << "Error! This point should be a center" << std::endl;

    auto entry = distance(voronoi_centers.begin(), it);
    vor.voronoi_regions[entry].push_back(i);
  }
}
vector<int> add_point_to_sampling(iVd &vor, const geodesic_solver &solver,
                                  const DrawableTrimesh<> &mesh, const int vid,
                                  const int type) {
  vor.distances[vid] = 0;
  vor.voronoi_tags[vid] = vid;

  vector<int> curr_voronoi = update_geodesic_distances_tracking_propoagation(
      vor.distances, solver, {vid}, type, vor.R);

  for (auto vert : curr_voronoi)
    vor.voronoi_tags[vert] = vid;

  return curr_voronoi;
}
iVd farthest_point_sampling(vector<int> &seeds, const geodesic_solver &solver,
                            const DrawableTrimesh<> &mesh, const int type,
                            const int k) {
  iVd vor = {};
  if (seeds.size() > 0) {
    vor.voronoi_tags = vector<int>(solver.graph.size(), seeds[0]);

    vor.region_normal_deviation.push_back({__DBL_MAX__, -1});

    vor.distances = compute_geodesic_distances(solver, {seeds[0]}, type);

    update_voronoi_regions(vor, seeds, mesh);
    vor.R = __DBL_MAX__;

    for (auto i = 1; i < seeds.size(); ++i) {

      auto curr_voronoi =
          add_point_to_sampling(vor, solver, mesh, seeds[i], type);
      update_voronoi_verts(vor, mesh, curr_voronoi);
      update_voronoi_regions(vor, seeds, mesh);
      vor.R = *std::max_element(vor.distances.begin(), vor.distances.end());
    }
  } else {
    seeds.reserve(k);
    seeds.push_back(0);
    vor.voronoi_tags = vector<int>(solver.graph.size(), 0);

    vor.region_normal_deviation.push_back({__DBL_MAX__, -1});

    vor.distances = compute_geodesic_distances(solver, {0}, type);

    update_voronoi_regions(vor, seeds, mesh);
    vor.R = __DBL_MAX__;

    for (auto i = 1; i < k; ++i) {

      int next_seed = std::distance(
          vor.distances.begin(),
          max_element(vor.distances.begin(), vor.distances.end()));
      auto curr_voronoi =
          add_point_to_sampling(vor, solver, mesh, next_seed, type);
      seeds.push_back(next_seed);
      update_voronoi_verts(vor, mesh, curr_voronoi);
      update_voronoi_regions(vor, seeds, mesh);
      vor.R = *std::max_element(vor.distances.begin(), vor.distances.end());
    }
  }
  return vor;
}