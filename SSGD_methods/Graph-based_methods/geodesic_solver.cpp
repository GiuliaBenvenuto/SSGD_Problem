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

  //solver.print_graph();
  return solver;
}


// ---------- Lanthier ----------
// uint add_node(geodesic_solver &solver) {
//     solver.graph.push_back(vector<geodesic_solver::graph_edge>());
//     return solver.graph.size() - 1; // Return the new node's index
// }


// void add_directed_arc(geodesic_solver &solver, uint na, uint nb, float w) {
//     if (na >= solver.graph.size() || nb >= solver.graph.size()) return;
//     geodesic_solver::graph_edge edge;
//     edge.node = nb;
//     edge.length = w; // Weight of the arc
//     solver.graph[na].push_back(edge);
// }


// void add_undirected_arc(geodesic_solver &solver, uint na, uint nb, float w) {
//     add_directed_arc(solver, na, nb, w);
//     add_directed_arc(solver, nb, na, w);
// }


// int add_node(geodesic_solver &solver, cinolib::vec3d p, uint deg) {
//     // solver.graph.push_back({});
//     // return solver.graph.size() - 1;
//     geodesic_solver::graph_edge edge;
//     edge.node = -1;  // Initialize to -1 since there are no adjacent nodes yet
//     edge.length = DBL_MAX;  // Initialize to maximum length
    
//     // Create a new adjacency list for the node
//     std::vector<geodesic_solver::graph_edge> adj_list(deg, edge);

//     // Add the new adjacency list to the graph
//     solver.graph.push_back(adj_list);
//     return solver.graph.size() - 1;
// }

// int add_node_between(geodesic_solver &solver, cinolib::vec3d p, uint deg, int vertex_index_1, int vertex_index_2) {
//     // Create a new node
//     geodesic_solver::graph_edge edge;
//     edge.node = -1;  // Initialize to -1 since there are no adjacent nodes yet
//     edge.length = DBL_MAX;  // Initialize to maximum length
    
//     // Create a new adjacency list for the node
//     std::vector<geodesic_solver::graph_edge> adj_list(deg, edge);

//     // Insert the new adjacency list at the appropriate position between the two vertices
//     int new_node_index = vertex_index_1 + 1; // Insert between the two vertices
//     solver.graph.insert(solver.graph.begin() + new_node_index, adj_list);

//     return new_node_index;
// }

// int add_node(geodesic_solver &solver, cinolib::vec3d p) {
//     // Create a new node
//     geodesic_solver::graph_edge edge;
//     edge.node = -1;  // Initialize to -1 since there are no adjacent nodes yet
//     edge.length = DBL_MAX;  // Initialize to maximum length
    
//     // Create a new adjacency list for the node
//     std::vector<geodesic_solver::graph_edge> adj_list;
//     adj_list.push_back(edge);

//     // Add the new adjacency list to the graph
//     solver.graph.push_back(adj_list);
//     return solver.graph.size() - 1;
// }

// int add_node_between(geodesic_solver &solver, cinolib::vec3d p, int vertex_index_1, int vertex_index_2) {
//     std::vector<geodesic_solver::graph_edge> adj_list; 
    

//     // Insert the new adjacency list at the appropriate position between the two vertices
//     int new_node_index = vertex_index_1 + 1; // Insert between the two vertices
//     solver.graph.insert(solver.graph.begin() + new_node_index, adj_list);

//     return new_node_index;
// }

// int add_node_prova(geodesic_solver &solver, cinolib::vec3d p, uint deg, size_t index) {
//     geodesic_solver::graph_edge edge;
//     edge.node = -1;  // No adjacent node yet
//     edge.length = DBL_MAX;

//     std::vector<geodesic_solver::graph_edge> adj_list;
//     adj_list.push_back(edge);

//     // Check if index is specified and within the current graph size or use push_back if not
//     if (index != std::string::npos && index < solver.graph.size()) {
//         solver.graph.reserve(deg);
//         solver.graph.insert(solver.graph.begin() + index, adj_list);
//     } else {
//         solver.graph.reserve(deg);
//         solver.graph.push_back(adj_list);
//     }
//     return index != std::string::npos ? index : solver.graph.size() - 1;
// }




void add_directed_arc(geodesic_solver &solver, int a, int b, float length) {
    if (a >= 0 && a < solver.graph.size() && b >= 0 && b < solver.graph.size()) {
      solver.graph[a].push_back({b, length});
    } else {
      cout << "ERROR: invalid node index" << endl;
    }
}


std::vector<cinolib::vec3d> sample_uniform(DrawableTrimesh<> & m, uint e, float step)
// sample points on edge e distributed uniformly along the edge:
{
    std::vector<cinolib::vec3d> samples;
    cinolib::vec3d v0 = m.edge_vert(e,0);
    cinolib::vec3d v1 = m.edge_vert(e,1);
    double d = v0.dist(v1);                 // edge length
    uint ns = (uint)(d/step+0.5);           // #steiner pts on this edge
    double mystep = d/ns;                   // actual step on this edge
    for (uint s=0;s<ns;s++)         // generate all Steiner points on e
    {
        samples.push_back( (1.0-(mystep+s*mystep)/d)*v0+((mystep+s*mystep)/d)*v1 );
    }
    return samples;
}


void sample_mesh_uniform(DrawableTrimesh<> & m, uint pxedge,
                         std::vector<cinolib::vec3d> & SteinerPoints, std::vector< std::pair< uint,uint > > & SteinerPerEdge)
{
    float step = m.edge_avg_length()/(pxedge+1);    // step on avg edge
    SteinerPoints.reserve(pxedge*m.num_edges());
    SteinerPerEdge.resize(m.num_edges());
    // sample Steiner points on edges
    for (uint e=0;e<m.num_edges();e++)
    {
        std::vector<cinolib::vec3d> e_samp = sample_uniform(m,e,step);  // generate all Steiner points on e
        uint fst = SteinerPoints.size();                                // index of first point on this edge
        for (uint s=0;s<e_samp.size();s++)                              // add all Steinter points
            SteinerPoints.push_back(e_samp[s]);
        SteinerPerEdge[e] = std::make_pair(fst,SteinerPoints.size()-fst);
    }

}


void sample_mesh_Steiner(DrawableTrimesh<> & m, uint pxedge,
                         std::vector<cinolib::vec3d> & SteinerPoints, std::vector< std::pair< uint,uint > > & SteinerPerEdge)
{ sample_mesh_uniform(m,pxedge,SteinerPoints,SteinerPerEdge); } 



void add_node(geodesic_solver &solver) {
  // geodesic_solver::graph_edge edge;
  // edge.node = -1;  // No adjacent node yet
  // edge.length = DBL_MAX;

  // std::vector<geodesic_solver::graph_edge> adj_list;
  // solver.graph.push_back(adj_list);

  std::vector<geodesic_solver::graph_edge> adj_list; 
  solver.graph.push_back(adj_list); 
  // print the node added

  //solver.graph.push_back({});
}


geodesic_solver compute_fine_graph(DrawableTrimesh<> &m, uint pxedge) {
    auto solver = geodesic_solver{};
    solver.graph.clear();

    std::vector<cinolib::vec3d> SteinerPoints;
    std::vector<std::pair<uint, uint>> SteinerPerEdge;
    sample_mesh_Steiner(m, pxedge, SteinerPoints, SteinerPerEdge);
    cout << "Steiner points: " << SteinerPoints.size() << endl;
    cout << "Steiner per edge: " << SteinerPerEdge.size() << endl;

    // Add nodes to graph: vertices + Steiner points
    uint offset = m.num_verts();
    cout << "offset: " << offset << endl;
    uint numn = offset + SteinerPoints.size();
    cout << "numn: " << numn << endl;
    solver.graph.reserve(numn);             // reserve number of nodes to speedup insertion

    for (uint i = 0; i < offset; i++) {
      add_node(solver);  // Add each vertex as a node
    }
    for (uint i = 0; i < SteinerPoints.size(); i++) {
      add_node(solver);
    }

    //solver.print_graph();

    
    //cout << "EDGES: " << m.num_edges() << endl;
    //cout << "VERTICES: " << offset << " STEINER POINTS: " << SteinerPoints.size() << endl;
    //cout << "Number of nodes: " << solver.graph.size() << endl;

    cout << "Number of polys: " << m.num_polys() << endl;
    
    for (uint i = 0; i < m.num_polys(); i++) {  // For all tris (face)
        for (uint j = 0; j < 3; j++) {          // For all edges of the face

            uint e = m.poly_edge_id(i, j);
            uint lastk = SteinerPerEdge[e].first + SteinerPerEdge[e].second;
            //cout << "lastk: " << lastk << endl;
            //cout << "SteinerPerEdge[e].first: " << SteinerPerEdge[e].first << endl;
            for (uint k = SteinerPerEdge[e].first; k < lastk; k++) {
                uint e1 = m.poly_edge_id(i, (j + 1) % 3);
                uint e2 = m.poly_edge_id(i, (j + 2) % 3);

                // Connect to Steiner points on the next edge
                uint lasth = SteinerPerEdge[e1].first + SteinerPerEdge[e1].second;
                //cout << "lasth: " << lasth << endl;
                for (uint h = SteinerPerEdge[e1].first; h < lasth; h++) {
                    float w = SteinerPoints[k].dist(SteinerPoints[h]);
                    // cout << "k: " << k << " h: " << h << " w: " << w << endl; 
                    add_directed_arc(solver, k + offset, h + offset, w);
                }

                // Connect to Steiner points on the following edge
                lasth = SteinerPerEdge[e2].first + SteinerPerEdge[e2].second;
                // cout << "lasth: " << lasth << endl;
                for (uint h = SteinerPerEdge[e2].first; h < lasth; h++) {
                    float w = SteinerPoints[k].dist(SteinerPoints[h]);
                    // cout << "k: " << k << " h: " << h << " w: " << w << endl; 
                    add_directed_arc(solver, k + offset, h + offset, w);
                }
            }
        }
    }
    //solver.print_graph();

    // Add arcs to vertices and along edges
    //cout << "Number of edges: " << m.num_edges() << endl;
    
    for (uint e = 0; e < m.num_edges(); e++) {
        // Get endpoints and opposite vertices
        uint v0, v1, v2, v3;
        v0 = m.edge_vert_id(e, 0);
        v1 = m.edge_vert_id(e, 1);
        //cout << "v0: " << v0 << " v1: " << v1 << endl;
        
        connect_nodes(solver, v0, v1, m.vert(v0).dist(m.vert(v1)));

        if (SteinerPerEdge[e].second == 0) continue;

        std::vector<uint> ff = m.adj_e2p(e);
        v2 = m.vert_opposite_to(ff[0], v0, v1);
        v3 = m.vert_opposite_to(ff[1], v0, v1);

        uint firstk = SteinerPerEdge[e].first;
        uint lastk = firstk + SteinerPerEdge[e].second - 1;

        float w0 = SteinerPoints[firstk].dist(m.vert(v0));
        float w1 = SteinerPoints[firstk].dist(m.vert(v1));
        //cout << "w0: " << w0 << " w1: " << w1 << endl;

        
        // Ensure v0 is closer to the first Steiner point
        if (w1 < w0) {
            int v=v0; v0=v1; v1=v;
            w0=w1;
        }

        // Connect the first and last Steiner points
        //cout << "firstk + offset: " << firstk + offset << endl;
        //cout << "lastk + offset: " << lastk + offset << endl;
        //cout << SteinerPoints[lastk].dist(m.vert(v1)) << endl;
        
        add_directed_arc(solver, firstk + offset, v0, w0);
        add_directed_arc(solver, lastk + offset, v1, SteinerPoints[lastk].dist(m.vert(v1)));

        // Connect first and last Steiner points to v2 and v3
        connect_nodes(solver, firstk + offset, v2, SteinerPoints[firstk].dist(m.vert(v2)));
        connect_nodes(solver, firstk + offset, v3, SteinerPoints[firstk].dist(m.vert(v3)));
        connect_nodes(solver, lastk + offset, v2, SteinerPoints[lastk].dist(m.vert(v2)));
        connect_nodes(solver, lastk + offset, v3, SteinerPoints[lastk].dist(m.vert(v3)));
        
        //cout << "SteinerPerEdge[e].second: " << SteinerPerEdge[e].second << endl;
        if (SteinerPerEdge[e].second > 1) {
            add_directed_arc(solver, firstk + offset, firstk + offset + 1, SteinerPoints[firstk].dist(SteinerPoints[firstk + 1]));
            add_directed_arc(solver, lastk + offset, lastk + offset - 1, SteinerPoints[lastk].dist(SteinerPoints[lastk - 1]));
        }
        firstk++;

        // Connect all other Steiner points
        //cout << "firstk: " << firstk << " lastk: " << lastk << endl;
        for (uint k = firstk; k < lastk; k++) {
            connect_nodes(solver, k + offset, v2, SteinerPoints[k].dist(m.vert(v2)));
            connect_nodes(solver, k + offset, v3, SteinerPoints[k].dist(m.vert(v3)));
            add_directed_arc(solver, k + offset, k + offset - 1, SteinerPoints[k].dist(SteinerPoints[k - 1]));
            add_directed_arc(solver, k + offset, k + offset + 1, SteinerPoints[k].dist(SteinerPoints[k + 1]));
        }
    }
    //solver.print_graph();
    return solver;

}


