#include "shortest_path.h"

int vert_offset(const DrawableTrimesh<> &m, const int tid, const int vid) {
  int k = -1;
  vector<uint> vids = m.poly_verts_id(tid);
  for (uint j = 0; j < 3; ++j)
    if (vids[j] == vid)
      return j;

  return k;
}
vec3d mesh_point_pos(const DrawableTrimesh<> &m, const mesh_point &p) {
  vector<vec3d> verts = m.poly_verts(p.tid);
  return interpolate(verts[0], verts[1], verts[2], p.bary);
}
vec2i oriented_edge(const DrawableTrimesh<> &m, const int tid,
                    const int next_tid) {
  uint eid = m.edge_shared(tid, next_tid);
  uint vid0 = m.edge_vert_id(eid, 0);
  uint vid1 = m.edge_vert_id(eid, 1);

  uint k = vert_offset(m, tid, vid0);
  if (m.poly_vert_id(tid, (k + 1) % 3) == vid1)
    return vec2i(vid0, vid1);
  else
    return vec2i(vid1, vid0);
}

void dump_path(const DrawableTrimesh<> &m, const vector<int> &strip,
               const vector<double> &lerps, const mesh_point &src,
               const mesh_point &tgt) {
  static int path_number = -1;
  cout << "dump current path in path_" + to_string(path_number + 1) +
              ".polylines.txt\n";
  ofstream out("path_" + to_string(++path_number) + ".polylines.txt");
  auto get_point_pos = [&m](const double lerp, const int curr_tid,
                            const int next_tid) {
    vec2i edge = oriented_edge(m, curr_tid, next_tid);
    return (1 - lerp) * m.vert(edge.x()) + lerp * m.vert(edge.y());
  };
  vec3d curr_pos = mesh_point_pos(m, src);
  out << strip.size() + 1 << " " << curr_pos.x() << " " << curr_pos.y() << " "
      << curr_pos.z();
  for (size_t i = 0; i < strip.size() - 1; ++i) {
    curr_pos = get_point_pos(lerps[i], strip[i], strip[i + 1]);
    out << " " << curr_pos.x() << " " << curr_pos.y() << " " << curr_pos.z();
  }
  curr_pos = mesh_point_pos(m, tgt);
  out << " " << curr_pos.x() << " " << curr_pos.y() << " " << curr_pos.z()
      << "\n";
}
void dump_portal(const vec2d &left, const vec2d &right) {

  ofstream out("portal.polylines.txt");

  out << 2 << " " << left.x() << " " << left.y() << " " << 0 << " " << right.x()
      << " " << right.y() << " " << 0 << "\n";
}
void dump_rigth_and_left_bound(const vec2d &left, const vec2d &right) {

  ofstream out("bound.polylines.txt");

  out << 2 << " " << left.x() << " " << left.y() << " " << 0 << " " << right.x()
      << " " << right.y() << " " << 0 << "\n";
}
void dump_points(const vector<funnel_point> &points) {

  ofstream out("points.pointset.txt");
  out << points.size();
  for (auto &point : points)
    out << " " << point.pos.x() << " " << point.pos.y() << " " << 0 << " ";
}
mesh_point get_point_from_vert(const DrawableTrimesh<> &m, const int vid) {
  uint tid = m.adj_v2p(vid)[0];
  int k = vert_offset(m, tid, vid);
  assert(k != -1);
  vec3d bary = vec3d(0., 0., 0.);
  bary[k] = 1.0;
  return mesh_point(tid, bary);
}
double cross(const vec2d &v0, const vec2d &v1) {
  return v0.x() * v1.y() - v0.y() * v1.x();
}
inline bool are_equal(const vec2d &v, const vec2d &w) {
  return (v.x() == w.x() && v.y() == w.y());
}
tuple<bool, int> point_is_edge(const mesh_point &p, const double &tol = 1e-5) {
  vec3d bary = p.bary;
  if (bary[0] > tol && bary[1] > tol && bary[2] <= tol)
    return {true, 0};
  if (bary[1] > tol && bary[2] > tol && bary[0] <= tol)
    return {true, 1};
  if (bary[2] > tol && bary[0] > tol && bary[1] <= tol)
    return {true, 2};

  return {false, -1};
}
bool check_strip(const DrawableTrimesh<> &m, const vector<int> &strip) {
  if (strip.size() == 0)
    return true;
  unordered_set<int> faces;
  faces.insert(strip[0]);
  for (size_t i = 1; i < strip.size(); ++i) {
    assert(faces.count(strip[i]) == 0);
    faces.insert(strip[i]);
    assert(m.polys_are_adjacent(strip[i - 1], strip[i]));
  }

  return true;
}

static tuple<bool, int> point_is_vert(const mesh_point &p,
                                      const double &tol = 1e-5) {
  vec3d bary = p.bary;

  if (bary[0] > tol && bary[1] <= tol && bary[2] <= tol)
    return {true, 0};
  if (bary[1] > tol && bary[0] <= tol && bary[2] <= tol)
    return {true, 1};
  if (bary[2] > tol && bary[0] <= tol && bary[1] <= tol)
    return {true, 2};

  return {false, -1};
}
tuple<bool, vec3d> point_in_triangle(const vec3d &point,
                                     const vector<vec3d> &verts,
                                     const double &tol = 1e-5) {
  // http://www.r-5.org/files/books/computers/algo-list/realtime-3d/Christer_Ericson-Real-Time_Collision_Detection-EN.pdf
  // pag.48
  vec3d b = vec3d(0., 0., 0.);
  vec3d v0 = verts[0];
  vec3d v1 = verts[1];
  vec3d v2 = verts[2];

  vec3d u = v1 - v0, v = v2 - v0, w = point - v0;
  double d00 = u.norm_sqrd(), d01 = u.dot(v), d11 = v.norm_sqrd(),
         d20 = w.dot(u), d21 = w.dot(v), d = d00 * d11 - d01 * d01;

  if (d == 0)
    return {false, b};

  b[2] = (d00 * d21 - d01 * d20) / d;
  assert(!isnan(b[2]));
  b[1] = (d11 * d20 - d01 * d21) / d;
  assert(!isnan(b[1]));
  b[0] = 1 - b[1] - b[2];
  assert(!isnan(b[0]));

  for (auto i = 0; i < 3; ++i) {
    if (b[i] < -tol || b[i] > 1.0 + tol)
      return {false, vec3d(0., 0., 0.)};
  }

  return {true, b};
}
double intersect_segments(const vec2d &start1, const vec2d &end1,
                          const vec2d &start2, const vec2d &end2) {
  if (are_equal(end1, start2))
    return 0;
  if (are_equal(end2, start1))
    return 1;
  if (are_equal(start2, start1)) {
    if (are_equal(end2, end1))
      return 1; // the portail and path segment coincide
    return 0;
  }
  if (are_equal(end2, end1))
    return 1;
  auto a = end1 - start1;   // direction of line a
  auto b = start2 - end2;   // direction of line b, reversed
  auto d = start2 - start1; // left-hand side
  auto det = a.x() * b.y() - a.y() * b.x();
  assert(det);
  return (a.x() * d.y() - a.y() * d.x()) / det;
}

vec2d intersect_circles(const vec2d &c2, const double &R2, const vec2d &c1,
                        const double &R1) {
  auto R = (c2 - c1).norm_sqrd();
  assert(R > 0);
  auto invR = double(1) / R;
  vec2d result = c1 + c2;

  result = result + (c2 - c1) * ((R1 - R2) * invR);
  auto A = 2 * (R1 + R2) * invR;
  auto B = (R1 - R2) * invR;
  auto s = A - B * B - 1;
  assert(s >= 0);
  result = result + vec2d(c2.y() - c1.y(), c1.x() - c2.x()) * sqrt(s);
  return result / 2.;
}

void clean_strip(vector<int> &strip, mesh_point &src, mesh_point &tgt,
                 const DrawableTrimesh<> &m) {

  if (strip.size() == 0)
    return;
  auto [source_is_vert, kvs] = point_is_vert(src);
  auto [source_is_edge, kes] = point_is_edge(src);
  auto [target_is_vert, kvt] = point_is_vert(tgt);
  auto [target_is_edge, ket] = point_is_edge(tgt);
  if (source_is_vert) {
    auto first_face = strip.begin();
    uint vid = m.poly_vert_id(src.tid, kvs);
    for (size_t i = 1; i < strip.size(); ++i) {
      int k = vert_offset(m, strip[i], vid);
      if (k == -1)
        break;
      ++first_face;
    }
    strip.erase(strip.begin(), first_face);
    if (strip.empty())
      return;

    int k = vert_offset(m, strip[0], vid);
    assert(k != -1);
    vec3d bary = vec3d(0., 0., 0.);
    bary[k] = 1;
    src = mesh_point(strip[0], bary);

  } else if (source_is_edge) {
    vec3d pos = mesh_point_pos(m, src);
    auto [is_in_next_tid, bary] =
        point_in_triangle(pos, m.poly_verts(strip[1]));
    if (is_in_next_tid) {
      strip.erase(strip.begin(), next(strip.begin()));
      src = mesh_point(strip[0], bary);
    }
  }

  if (target_is_vert) {
    auto last_face = strip.end();
    uint vid = m.poly_vert_id(tgt.tid, kvt);
    for (auto rit = strip.rbegin(); rit != strip.rend(); ++rit) {
      if (rit == strip.rend() - 1) {
        kvt = vert_offset(m, *rit, vid);
        assert(kvt != -1);
        strip.clear();
        return;
      }
      kvt = vert_offset(m, *(rit + 1), vid);
      if (kvt == -1)
        break;
      --last_face;
    }

    strip.erase(last_face, strip.end());

    if (strip.empty())
      return;

    int k = vert_offset(m, strip.back(), vid);
    vec3d bary = vec3d(0., 0., 0.);
    bary[k] = 1;
    tgt = mesh_point(strip.back(), bary);
  } else if (target_is_edge) {
    vec3d pos = mesh_point_pos(m, tgt);
    auto [is_in_prev_tid, bary] =
        point_in_triangle(pos, m.poly_verts(strip.rbegin()[1]));
    if (is_in_prev_tid) {
      strip.pop_back();
      tgt = mesh_point(strip.back(), bary);
    }
  }
}
array<vec2d, 2> init_source_triangle(const DrawableTrimesh<> &m,
                                     const mesh_point &src, const uint eid) {

  uint vid0 = m.edge_vert_id(eid, 0);
  uint vid1 = m.edge_vert_id(eid, 1);
  int k = vert_offset(m, src.tid, vid0);
  assert(k != -1);
  if (m.poly_vert_id(src.tid, (k + 1) % 3) != vid1) {
    k = (k + 2) % 3;
    swap(vid0, vid1);
  }
  uint vid2 = m.vert_opposite_to(src.tid, vid0, vid1);

  vec2d flat_v0 = vec2d(0., 0.);
  vec2d flat_v1 = vec2d(0, (m.vert(vid1) - m.vert(vid0)).norm());
  double r0 = (m.vert(vid2) - m.vert(vid0)).norm_sqrd();
  double r1 = (m.vert(vid2) - m.vert(vid1)).norm_sqrd();
  vec2d flat_v2 = intersect_circles(flat_v0, r0, flat_v1, r1);
  vec2d point_coords = flat_v0 * src.bary[k] + flat_v1 * src.bary[(k + 1) % 3] +
                       flat_v2 * src.bary[(k + 2) % 3];
  array<vec2d, 2> res;
  res[0] = flat_v0 - point_coords;
  res[1] = flat_v1 - point_coords;
  return res;
}
array<vec2d, 2> init_target_triangle(const DrawableTrimesh<> &m,
                                     const array<vec2d, 2> &prev,
                                     const mesh_point &tgt,
                                     const int prev_tid) {
  vec2i edge = oriented_edge(m, prev_tid, tgt.tid);

  int vid0 = edge.x();
  int vid1 = edge.y();

  int k = vert_offset(m, tgt.tid, vid0);
  uint v = m.vert_opposite_to(tgt.tid, vid0, vid1);
  double r0 = (m.vert(vid0) - m.vert(v)).norm_sqrd();
  double r1 = (m.vert(vid1) - m.vert(v)).norm_sqrd();
  vec2d flat_v = intersect_circles(prev[1], r1, prev[0], r0);
  vec2d tgt_coords = prev[0] * tgt.bary[k] + flat_v * tgt.bary[(k + 1) % 3] +
                     prev[1] * tgt.bary[(k + 2) % 3];

  array<vec2d, 2> res;
  res[0] = tgt_coords;
  res[1] = tgt_coords;

  return res;
}
array<vec2d, 2> unfold_face(const DrawableTrimesh<> &m, const int curr_tid,
                            const int next_tid, const int next_eid,
                            const array<vec2d, 2> &flat_tid) {
  vec2i edge = oriented_edge(m, curr_tid, next_tid);
  int vid0 = edge.x();
  int vid1 = edge.y();
  uint v = m.vert_opposite_to(next_tid, vid0, vid1);
  double r0 = (m.vert(v) - m.vert(vid0)).norm_sqrd();
  double r1 = (m.vert(v) - m.vert(vid1)).norm_sqrd();
  vec2d v_flat = intersect_circles(flat_tid[1], r1, flat_tid[0], r0);

  uint right = m.edge_id(vid1, v);
  uint left = m.edge_id(vid0, v);
  array<vec2d, 2> res;
  if (next_eid == right) {
    res[0] = v_flat;
    res[1] = flat_tid[1];
  } else if (next_eid == left) {
    res[0] = flat_tid[0];
    res[1] = v_flat;
  } else
    assert(false);

  return res;
}
vector<array<vec2d, 2>> unfold_strip(const DrawableTrimesh<> &m,
                                     const vector<int> &strip,
                                     const mesh_point &src,
                                     const mesh_point &tgt) {
  assert(strip[0] == src.tid && strip.back() == tgt.tid);
  vector<array<vec2d, 2>> result(strip.size());

  result[0] = init_source_triangle(m, src, m.edge_shared(strip[0], strip[1]));
  for (size_t i = 1; i < strip.size() - 1; ++i)
    result[i] =
        unfold_face(m, strip[i - 1], strip[i],
                    m.edge_shared(strip[i], strip[i + 1]), result[i - 1]);

  result.back() =
      init_target_triangle(m, result.rbegin()[1], tgt, strip.rbegin()[1]);

  return result;
}
int max_curvature_point(const vector<funnel_point> &path) {
  // Among vertices around which the path curves, find the vertex
  // with maximum angle. We are going to fix that vertex. Actually, max_index is
  // the index of the first face containing that vertex.
  size_t max_index = -1;
  double max_angle = 0.;
  for (size_t i = 1; i < path.size() - 1; ++i) {
    vec2d pos = path[i].pos;
    vec2d prev = path[i - 1].pos;
    vec2d next = path[i + 1].pos;

    vec2d v0 = pos - prev;
    vec2d v1 = next - pos;

    double angle = v0.angle_rad(v1);
    if (angle > max_angle) {
      max_index = path[i].face;
      max_angle = angle;
    }
  }

  return max_index;
}
vector<double> funnel(const vector<array<vec2d, 2>> &portals,
                      size_t &max_index) {
  // Find straight path.
  vec2d start = vec2d(0., 0.);
  int apex_index = 0;
  int rigth_index = 0;
  int left_index = 0;
  vec2d apex = start;
  vec2d right_bound = portals[0][0];
  vec2d left_bound = portals[0][1];

  // Add start point.
  vector<funnel_point> points = vector<funnel_point>{{apex_index, apex}};
  points.reserve(portals.size());

  auto area = [](const vec2d a, const vec2d b, const vec2d c) {
    return cross(b - a, c - a);
  };

  for (size_t i = 0; i < portals.size(); ++i) {
    auto right = portals[i][0], left = portals[i][1];

    // Update left vertex.
    if (area(apex, left_bound, left) <= 0) {
      if (apex == left_bound || area(apex, right_bound, left) > 0) {
        // Tighten the funnel.
        left_bound = left;
        left_index = i;

      } else {
        // Right over right, insert right to path and restart scan from
        // portal right point.
        if (!are_equal(right_bound, apex)) {
          points.push_back({rigth_index, right_bound});
          // Make current right the new apex.
          apex = right_bound;
          apex_index = rigth_index;
          // Reset portal
          right_bound = apex;
          left_bound = apex;
          rigth_index = apex_index;
          left_index = apex_index;
          // Restart scan
          i = apex_index;

          continue;
        }
      }
    }

    // Update right vertex.
    if (area(apex, right_bound, right) >= 0) {
      if (apex == right_bound || area(apex, left_bound, right) < 0) {
        // Tighten the funnel.
        right_bound = right;
        rigth_index = i;

      } else {
        if (!are_equal(left_bound, apex)) {
          points.push_back({left_index, left_bound});
          // Make current left the new apex.
          apex = left_bound;
          apex_index = left_index;
          // Reset portal
          right_bound = apex;
          left_bound = apex;
          rigth_index = apex_index;
          left_index = apex_index;
          // Restart scan
          i = apex_index;

          continue;
        }
      }
    }
  }

  // This happens when we got an apex on the last edge of the strip
  if (!are_equal(points.back().pos, portals.back()[0])) {
    points.push_back({(int)portals.size() - 1, portals.back()[0]});
  }
  assert(points.back().pos == portals.back()[0]);
  assert(points.back().pos == portals.back()[1]);

  vector<double> lerps;
  lerps.reserve(portals.size());
  for (size_t i = 0; i < points.size() - 1; i++) {
    auto a = points[i].pos;
    auto b = points[i + 1].pos;
    for (auto k = points[i].face; k < points[i + 1].face; ++k) {
      auto portal = portals[k];
      double s = intersect_segments(a, b, portal[0], portal[1]);
      lerps.push_back(std::clamp(s, 0.0, 1.0));
    }
  }

  auto index = 1;

  for (size_t i = 0; i < portals.size(); ++i) {
    if ((portals[i][0] == points[index].pos) ||
        (portals[i][1] == points[index].pos)) {
      points[index].face = i;
      index += 1;
    }
  }

  max_index = max_curvature_point(points);
  assert(lerps.size() == portals.size() - 1);
  return lerps;
}

vector<int> triangle_fan(const DrawableTrimesh<> &m, const int start_face,
                         const int end_face, const int vid, const bool CCW,
                         const vector<int> &strip, size_t &index) {
  vector<int> result;
  result.reserve(m.adj_v2p(vid).size());
  result.push_back(start_face);
  uint k = vert_offset(m, start_face, vid);
  uint vid0 = (CCW) ? m.poly_vert_id(start_face, (k + 2) % 3)
                    : m.poly_vert_id(start_face, (k + 1) % 3);
  uint eid = m.edge_id(vid, vid0);
  uint next_face = m.polys_adjacent_along(start_face, eid)[0];

  while (index >= 1 && strip[index - 1] == next_face) {
    result.pop_back();
    result.push_back(next_face);
    vid0 = m.vert_opposite_to(next_face, vid, vid0);
    eid = m.edge_id(vid, vid0);
    next_face = m.polys_adjacent_along(next_face, eid)[0];
    --index;
  }

  result.push_back(next_face);
  while (next_face != end_face) {

    vid0 = m.vert_opposite_to(next_face, vid, vid0);
    eid = m.edge_id(vid, vid0);
    next_face = m.polys_adjacent_along(next_face, eid)[0];
    result.push_back(next_face);
  }

  return result;
}
void straighten_path(vector<array<vec2d, 2>> &portals, vector<double> &lerps,
                     vector<int> &strip, mesh_point &src, mesh_point &tgt,
                     size_t index, const DrawableTrimesh<> &m) {
  int vertex = -1;

  // while(true) { this may never break...
  for (auto i = 0; i < strip.size() * 2 && index != -1; i++) {
    int new_vertex = -1;
    int face = strip[index];
    int next = strip[index + 1];
    vec2i edge = oriented_edge(m, face, next);
    bool CCW = false;
    if (lerps[index] == 0) {
      new_vertex = edge.x();
      CCW = true;
    } else if (lerps[index] == 1) {
      new_vertex = edge.y();
    }
    if (new_vertex == vertex)
      break;
    vertex = new_vertex;

    int curr_index = index + 1;

    int target_face;
    if (curr_index == strip.size())
      target_face = tgt.tid;
    else {
      while (vert_offset(m, strip[curr_index], vertex) != -1) {
        if (curr_index == strip.size() - 1) {
          target_face = tgt.tid;
          curr_index = strip.size();
          break;
        }

        ++curr_index;
      }

      target_face = strip[curr_index - 1];
    }

    vector<int> new_faces =
        triangle_fan(m, face, target_face, vertex, CCW, strip, index);
    vector<int> tmp(strip.begin(), strip.begin() + index);
    tmp.insert(tmp.end(), new_faces.begin(), new_faces.end());
    tmp.insert(tmp.end(), strip.begin() + curr_index, strip.end());
    strip.swap(tmp);
    clean_strip(strip, src, tgt, m);
    if (strip.size() == 0)
      return;
    portals = unfold_strip(m, strip, src, tgt);
    lerps = funnel(portals, index);
    assert(lerps.size() == strip.size() - 1);
  }
}

vector<int> k_ring(const DrawableTrimesh<> &m, const int vid, const int k,
                   const bool boundary_only) {

  vector<int> active_set = {vid};
  vector<int> ring;

  if (!boundary_only)
    ring.push_back(vid);
  for (uint i = 0; i < k; ++i) {
    vector<int> next_active_set;
    for (uint curr : active_set)
      for (uint nbr : m.vert_ordered_verts_link(curr)) {
        if (find(ring.begin(), ring.end(), nbr) == ring.end() && nbr != vid) {
          next_active_set.push_back(nbr);
          ring.push_back(nbr);
        }
      }

    active_set = next_active_set;
  }

  return ring;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// Geodesic Solver
array<vec2d, 3> init_flat_triangle(const DrawableTrimesh<> &m, const uint tid) {

  array<vec2d, 3> result;
  vector<vec3d> verts = m.poly_verts(tid);
  result[0] = vec2d(0., 0.);
  result[1] = vec2d(0., (verts[1] - verts[0]).norm());
  double r0 = (verts[2] - verts[0]).norm_sqrd();
  double r1 = (verts[2] - verts[1]).norm_sqrd();
  result[2] = intersect_circles(result[0], r0, result[1], r1);

  return result;
}
array<vec2d, 3> unfold_face(const DrawableTrimesh<> &m, const uint tid,
                            const uint adj, const array<vec2d, 3> &flat_tid) {

  array<vec2d, 3> result;
  uint eid = m.edge_shared(tid, adj);
  uint vid0 = m.edge_vert_id(eid, 0);
  uint vid1 = m.edge_vert_id(eid, 1);
  int k = vert_offset(m, tid, vid0);

  if (k == -1)
    std::cout << "You are messing something up" << std::endl;

  if (m.poly_vert_id(tid, (k + 1) % 3) != vid1) {
    k = (k + 2) % 3;
    swap(vid0, vid1);
  }
  uint opp = m.vert_opposite_to(adj, vid0, vid1);
  double r1 = (m.vert(vid1) - m.vert(opp)).norm_sqrd();
  double r0 = (m.vert(vid0) - m.vert(opp)).norm_sqrd();

  vec2d flat_opp =
      intersect_circles(flat_tid[(k + 1) % 3], r1, flat_tid[k], r0);

  int h = vert_offset(m, adj, vid0);
  result[h] = flat_tid[k];
  result[(h + 1) % 3] = flat_opp;
  result[(h + 2) % 3] = flat_tid[(k + 1) % 3];

  return result;
}

dual_geodesic_solver make_dual_geodesic_solver(const DrawableTrimesh<> &m) {
  uint F = m.num_polys();
  dual_geodesic_solver result;
  result.graph.resize(F);
  auto compute_dual_weights = [&m](const int tid0, const int tid1) {
    array<vec2d, 3> flat_tid = init_flat_triangle(m, tid0);
    array<vec2d, 3> flat_nei = unfold_face(m, tid0, tid1, flat_tid);

    vec2d c0 = interpolate(flat_tid, vec3d(1.0 / 3, 1.0 / 3, 1.0 / 3));
    vec2d c1 = interpolate(flat_nei, vec3d(1.0 / 3, 1.0 / 3, 1.0 / 3));

    return (c1 - c0).norm();
  };

  for (uint i = 0; i < F; ++i) {
    for (uint k = 0; k < 3; ++k) {
      uint adj = m.adj_p2p(i)[k];
      result.graph[i][k].node = adj;
      if (adj == -1)
        result.graph[i][k].length = DBL_MAX;
      else
        result.graph[i][k].length = compute_dual_weights(i, adj);
    }
  }

  return result;
}
template <typename Update, typename Stop, typename Exit>
void search_strip(vector<double> &field, vector<bool> &in_queue,
                  const dual_geodesic_solver &solver,
                  const DrawableTrimesh<> &m, int start, int end,
                  Update &&update, Stop &&stop, Exit &&exit) {
  auto destination_pos =
      interpolate(m.poly_verts(end), vec3d{1.0 / 3, 1.0 / 3, 1.0 / 3});

  auto estimate_dist = [&](int face) {
    vec3d p = interpolate(m.poly_verts(face), vec3d{1.0 / 3, 1.0 / 3, 1.0 / 3});
    return (p - destination_pos).norm();
  };
  field[start] = estimate_dist(start);

  // Cumulative weights of elements in queue. Used to keep track of the
  // average weight of the queue.
  double cumulative_weight = 0.0;

  // setup queue
  auto queue = std::deque<int>{};
  in_queue[start] = true;
  cumulative_weight += field[start];
  queue.push_back(start);

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
      auto neighbor = solver.graph[node][i].node;
      if (neighbor == -1)
        continue;

      // Distance of neighbor through this node
      auto new_distance = field[node];
      new_distance += solver.graph[node][i].length;
      new_distance += estimate_dist(neighbor);
      new_distance -= estimate_dist(node);

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
      if (update(node, neighbor, new_distance))
        return;
    }
  }
}
vector<int> compute_strip_tlv(const DrawableTrimesh<> &m,
                              const dual_geodesic_solver &solver, int start,
                              int end) {
  if (start == end)
    return {start};

  thread_local static auto parents = vector<int>{};
  thread_local static auto field = vector<double>{};
  thread_local static auto in_queue = vector<bool>{};

  if (parents.size() != solver.graph.size()) {
    parents.assign(solver.graph.size(), -1);
    field.assign(solver.graph.size(), DBL_MAX);
    in_queue.assign(solver.graph.size(), false);
  }

  // initialize once for all and sparsely cleanup at the end of every solve
  auto visited = vector<int>{start};
  auto sources = vector<int>{start};
  auto update = [&visited, end](int node, int neighbor, float new_distance) {
    parents[neighbor] = node;
    visited.push_back(neighbor);
    return neighbor == end;
  };
  auto stop = [](int node) { return false; };
  auto exit = [](int node) { return false; };

  search_strip(field, in_queue, solver, m, start, end, update, stop, exit);

  // extract_strip
  auto strip = vector<int>{};
  auto node = end;
  strip.reserve((int)sqrt(parents.size()));
  while (node != -1) {
    strip.push_back(node);
    node = parents[node];
  }

  // cleanup buffers
  for (auto &v : visited) {
    parents[v] = -1;
    field[v] = DBL_MAX;
    in_queue[v] = false;
  }
  // assert(check_strip(geometry.adjacencies, strip));
  return strip;
}
vector<vec3d> shortest_path(mesh_point &src, mesh_point &tgt,
                            const DrawableTrimesh<> &m,
                            const dual_geodesic_solver &solver) {
  auto get_point_pos = [&m](const double lerp, const int curr_tid,
                            const int next_tid) {
    vec2i edge = oriented_edge(m, curr_tid, next_tid);
    return (1 - lerp) * m.vert(edge.x()) + lerp * m.vert(edge.y());
  };
  vector<int> strip = compute_strip_tlv(m, solver, tgt.tid, src.tid);
  // strip_on_dual_graph(solver, tgt.tid, src.tid);
  //  reversed so the strip goes from src to tgt
  assert(check_strip(m, strip));
  clean_strip(strip, src, tgt, m);
  if (strip.size() < 2)
    return {mesh_point_pos(m, src), mesh_point_pos(m, tgt)};

  vector<array<vec2d, 2>> portals = unfold_strip(m, strip, src, tgt);
  size_t max_index = 0;
  vector<double> lerps = funnel(portals, max_index);

  straighten_path(portals, lerps, strip, src, tgt, max_index, m);

  vector<vec3d> result(strip.size() + 1);
  result[0] = mesh_point_pos(m, src);
  for (uint i = 1; i < strip.size(); ++i)
    result[i] = get_point_pos(lerps[i - 1], strip[i - 1], strip[i]);

  result.back() = mesh_point_pos(m, tgt);

  return result;
}

double path_length(const vector<vec3d> &path) {
  double len = 0;
  for (size_t i = 0; i < path.size() - 1; ++i)
    len += (path[i + 1] - path[i]).norm();

  return len;
}
