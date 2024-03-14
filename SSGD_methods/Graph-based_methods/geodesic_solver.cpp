#include "geodesic_solver.h"
#include <deque>
template <typename Update, typename Stop, typename Exit>
void visit_geodesic_graph(vector<double> &field, const geodesic_solver &solver,
                          const vector<int> &sources, Update &&update,
                          Stop &&stop, Exit &&exit) {
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
      double new_distance = field[node] + solver.graph[node][i].length;

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
    }
  }
}
template <typename Update, typename Stop, typename Exit>
void visit_geodesic_graph(vector<double> &field,
                          const dual_geodesic_solver &solver,
                          const vector<int> &sources, Update &&update,
                          Stop &&stop, Exit &&exit) {
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
      double new_distance = field[node] + solver.graph[node][i].length;

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
      update(node, neighbor);
    }
  }
}

// update the distance field "distances" by adding "sources"
// as futher sources and propagating with Dijkstra
void update_geodesic_distances(vector<double> &distances,
                               const geodesic_solver &solver,
                               const vector<int> &sources,
                               double max_distance) {

  auto update = [](int node, int neighbor) {};
  auto stop = [&](int node) { return distances[node] > max_distance; };
  auto exit = [](int node) { return false; };
  for (auto source : sources)
    distances[source] = 0.0;
  visit_geodesic_graph(distances, solver, sources, update, stop, exit);
}

// compute the distance field from "sources" with Dijkstra
vector<double> compute_geodesic_distances(const geodesic_solver &solver,
                                          const vector<int> &sources) {

  auto field = vector<double>(solver.graph.size(), DBL_MAX);
  for (auto source : sources)
    field[source] = 0.0;

  update_geodesic_distances(field, solver, sources);

  return field;
}
vector<int> strip_on_dual_graph(const dual_geodesic_solver &solver,
                                const int start, const int end) {
  if (start == end)
    return {};
  vector<int> parents(solver.graph.size(), -1);
  vector<double> field(solver.graph.size(), DBL_MAX);
  vector<int> sources = {start};
  auto update = [&parents, end](int node, int neighbor) {
    parents[neighbor] = node;
  };

  auto stop = [](int node) { return false; };
  auto exit = [&end](int node) { return node == end; };
  field[start] = 0;
  visit_geodesic_graph(field, solver, {start}, update, stop, exit);

  vector<int> strip;
  strip.reserve((int)sqrt(parents.size()));
  int node = end;
  while (node != -1) {
    strip.push_back(node);
    node = parents[node];
  }
  return strip;
}
//:::::::::::::::::::::::::::::::::::::: Primal solvers
void connect_nodes(geodesic_solver &solver, int a, int b, float length) {
  solver.graph[a].push_back({b, length});
  solver.graph[b].push_back({a, length});
}

double opposite_nodes_arc_length(const vector<vec3d> &positions, int a, int c,
                                 const vec2i &edge) {
  // Triangles (a, b, d) and (b, d, c) are connected by (b, d) edge
  // Nodes a and c must be connected.

  auto b = edge.x(), d = edge.y();
  auto ba = positions[a] - positions[b];
  auto bc = positions[c] - positions[b];
  auto bd = positions[d] - positions[b];
  ba.normalize();
  bd.normalize();
  auto cos_alpha = ba.dot(bd);
  auto cos_beta = bc.dot(bd);
  auto sin_alpha = sqrt(max(0.0, 1 - cos_alpha * cos_alpha));
  auto sin_beta = sqrt(max(0.0, 1 - cos_beta * cos_beta));

  // cos(alpha + beta)
  auto cos_alpha_beta = cos_alpha * cos_beta - sin_alpha * sin_beta;
  if (cos_alpha_beta <= -1)
    return DBL_MAX;

  // law of cosines (generalized Pythagorean theorem)
  ba = positions[a] - positions[b];
  bc = positions[c] - positions[b];
  bd = positions[d] - positions[b];
  auto len =
      ba.dot(ba) + bc.dot(bc) - ba.norm() * bc.norm() * 2 * cos_alpha_beta;

  if (len <= 0)
    return DBL_MAX;
  else
    return sqrt(len);
}

static void connect_opposite_nodes(geodesic_solver &solver,
                                   const vector<vec3d> &positions,
                                   const vector<uint> &tr0,
                                   const vector<uint> &tr1, const vec2i &edge) {
  auto opposite_vertex = [](const vector<uint> &tr, const vec2i &edge) -> int {
    for (auto i = 0; i < 3; ++i) {
      if (tr[i] != edge.x() && tr[i] != edge.y())
        return tr[i];
    }
    return -1;
  };

  auto v0 = opposite_vertex(tr0, edge);
  auto v1 = opposite_vertex(tr1, edge);
  if (v0 == -1 || v1 == -1)
    return;
  auto length = opposite_nodes_arc_length(positions, v0, v1, edge);
  connect_nodes(solver, v0, v1, length);
}

geodesic_solver make_geodesic_solver(const DrawableTrimesh<> &m,
                                     const bool geo_tangle) {
  auto solver = geodesic_solver{};
  solver.graph.resize(m.num_verts());

  int edge_count = 0;

  for (auto face = 0; face < m.num_polys(); face++) {
    for (auto k = 0; k < 3; k++) {
      auto a = m.poly_vert_id(face, k);
      auto b = m.poly_vert_id(face, (k + 1) % 3);

      // connect mesh edges
      auto len = (m.vert(a) - m.vert(b)).norm();
      if (a < b)
        connect_nodes(solver, a, b, len);
      edge_count++;

      if (!geo_tangle)
        continue;
      // connect opposite nodes
      auto neighbor = m.adj_p2p(face)[k];
      if (face < neighbor) {
        connect_opposite_nodes(solver, m.vector_verts(), m.adj_p2v(face),
                               m.adj_p2v(neighbor), vec2i{(int)a, (int)b});
        edge_count++;
      }
    }
  }
  if (geo_tangle)
    std::cout << "Number of edges GEOTANGLE: " << edge_count << std::endl;
  else
    std::cout << "Number of edges edge graph: " << edge_count << std::endl;
  return solver;
}