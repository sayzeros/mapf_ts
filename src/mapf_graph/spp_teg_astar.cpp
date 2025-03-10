#include "spp_teg.h"
#include "sy_log.h"

#include <queue>

#define SPP_DIST_BIG_VALUE 1e9

int b_horizon_margin = 5;


struct AstarLabel
{
public:
  shared_ptr<Vertex> vertex;
  int time;
  double g; // true distance from start
  double f; // approximated total distance from start to vertex to terminal, g+h
  bool reserved;
  bool reserved2;

  AstarLabel(shared_ptr<Vertex> &vertex,
             int time,
             double g,
             double f) : vertex(vertex), time(time), g(g), f(f), reserved(true), reserved2(true) {}

  AstarLabel(shared_ptr<Vertex> &vertex,
             int time,
             double g,
             double f,
             bool reserved) : vertex(vertex), time(time), g(g), f(f), reserved(reserved), reserved2(true) {}

  AstarLabel(shared_ptr<Vertex> &vertex,
             int time,
             double g,
             double f,
             bool reserved,
             bool reserved2) : vertex(vertex), time(time), g(g), f(f), reserved(reserved), reserved2(reserved2) {}
};

struct AstarLabelCompare
{
  bool operator()(const AstarLabel &a, const AstarLabel &b) const
  {
    return (a.f > b.f) ||
           (a.f == b.f && a.reserved > b.reserved) ||
           (a.f == b.f && a.reserved == b.reserved && a.reserved2 > b.reserved2) ||
           (a.f == b.f && a.reserved == b.reserved && a.reserved2 == b.reserved2 && a.g < b.g) ||
           (a.f == b.f && a.reserved == b.reserved && a.reserved2 == b.reserved2 && a.g == b.g);
  }
};

// Add mechanism to remove duplicated labels of a vertex to priority_queue
// if MANUALLY_REMOVE_DUPLICATED_LABEL, check it when pushed (unimplemented in project MAPF-TS)
// Else, check it when poped
class AstarQueue : public std::priority_queue<AstarLabel, vector<AstarLabel>, AstarLabelCompare>
{
public:
#ifdef MANUALLY_REMOVE_DUPLICATED_LABEL
  void push(AstarLabel &&label)
  {
    if (vertex_is_in_queue[label.vertex])
    {
      auto it = find(label.vertex);
      if (!AstarLabelCompare()(*it, label))
      {
        return;
      }
      remove(it);
    }
    else
    {
      vertex_is_in_queue[label.vertex] = true;
    }
    std::priority_queue<AstarLabel, vector<AstarLabel>, AstarLabelCompare>::push(label);
  }
#endif

  void pop()
  {
    const AstarLabel &top = std::priority_queue<AstarLabel, vector<AstarLabel>, AstarLabelCompare>::top();
#ifdef MANUALLY_REMOVE_DUPLICATED_LABEL
    vertex_is_in_queue[top.vertex] = false;
#endif
    vt_poped[make_tuple(top.vertex, top.time)] = true;

    std::priority_queue<AstarLabel, vector<AstarLabel>, AstarLabelCompare>::pop();
  }

#ifndef MANUALLY_REMOVE_DUPLICATED_LABEL
  AstarLabel top()
  {
    const AstarLabel &top = std::priority_queue<AstarLabel, vector<AstarLabel>, AstarLabelCompare>::top();

    if (vt_poped.contains(make_tuple(top.vertex, top.time)))
    {
      if (std::priority_queue<AstarLabel, vector<AstarLabel>, AstarLabelCompare>::empty())
        return top;
      std::priority_queue<AstarLabel, vector<AstarLabel>, AstarLabelCompare>::pop();
      return AstarQueue::top();
    }

    return top;
  }
#endif

  bool is_poped(AstarLabel &label)
  {
    return vt_poped.contains(make_tuple(label.vertex, label.time));
  }
  bool is_poped(shared_ptr<Vertex> vertex, int time)
  {
    return vt_poped.contains(make_tuple(vertex, time));
  }

private:
  sy_map<tuple<shared_ptr<Vertex>, int>, bool> vt_poped;
#ifdef MANUALLY_REMOVE_DUPLICATED_LABEL
  sy_map<shared_ptr<LSPLayeredVertex>, bool> vertex_is_in_queue;

  // return iterator with same vertex with input
  vector<AstarLabel>::iterator find(shared_ptr<LSPLayeredVertex> vertex)
  {
    for (auto it = this->c.begin(); it != this->c.end(); it += 1)
    {
      if (it->vertex == vertex)
      {
        return it;
      }
    }
    return this->c.end();
  }

  bool remove(vector<AstarLabel>::iterator it)
  {
    if (it == this->c.end())
    {
      return false;
    }
    if (it == this->c.begin())
    {
      std::priority_queue<AstarLabel, vector<AstarLabel>, AstarLabelCompare>::pop();
    }
    else
    {
      // remove element and re-heap
      this->c.erase(it);
      std::make_heap(this->c.begin(), this->c.end(), this->comp);
    }
    return true;
  }
#endif
};

/////////////////////////
/////////////////////////
/////////////////////////
/////////////////////////
/////////////////////////
/////////////////////////
/////////////////////////
/////////////////////////
/////////////////////////
/////////////////////////
/////////////////////////
/////////////////////////
/////////////////////////
/////////////////////////
/////////////////////////
/////////////////////////
/////////////////////////
/////////////////////////

sy_map<tuple<shared_ptr<Vertex>, shared_ptr<Vertex>>, double> plain_distance;

struct vd
{
  shared_ptr<Vertex> v;
  double dist;

  vd(shared_ptr<Vertex> v, double dist) : v(v), dist(dist) {}
};

struct cmp_vd
{
  bool operator()(vd l, vd r) { return l.dist > r.dist; }
};

double solve_spp_dijkstra(Graph &graph, shared_ptr<Vertex> vertex, shared_ptr<Vertex> goal)
{
  // sy_map<tuple<int, int>, shared_ptr<Arc>> pred;        // key: (vertex_id, time)
  sy_map<tuple<int, int>, double> dist_minus_big_value; // key: (vertex_id, time)

  std::priority_queue<vd, vector<vd>, cmp_vd> queue;
  sy_map<tuple<int, int>, bool> vt_visited; // default is false

  /* Dijkstra search */
  dist_minus_big_value[make_tuple(goal->id, 0)] = -SPP_DIST_BIG_VALUE;
  queue.push(vd(goal, dist_minus_big_value[make_tuple(goal->id, 0)]));

  while (!queue.empty())
  {
    auto top = queue.top();
    queue.pop();
    auto u = top.v;
    if (vt_visited[make_tuple(u->id, 0)])
      continue;
    vt_visited[make_tuple(u->id, 0)] = true;

    plain_distance[make_tuple(u, goal)] = dist_minus_big_value[make_tuple(u->id, 0)] + SPP_DIST_BIG_VALUE;

    if (u == vertex &&
        dist_minus_big_value[make_tuple(u->id, 0)] != 0)
    {
      break;
    }

    if (dist_minus_big_value[make_tuple(u->id, 0)] != 0)
    {
      for (auto &arc : graph.get_inin_arcs(u))
      {
        auto &v = arc->from;
        if (dist_minus_big_value[make_tuple(u->id, 0)] + arc->cost < dist_minus_big_value[make_tuple(v->id, 0)])
        {
          dist_minus_big_value[make_tuple(v->id, 0)] = dist_minus_big_value[make_tuple(u->id, 0)] + arc->cost;
          queue.push(vd(v, dist_minus_big_value[make_tuple(v->id, 0)]));
        }
      }
    }
  }

  return dist_minus_big_value[make_tuple(vertex->id, 0)] + SPP_DIST_BIG_VALUE;
}

double solve_spp_dijkstra_b(Graph &graph, shared_ptr<Vertex> vertex, shared_ptr<Vertex> goal)
{
  // sy_map<tuple<int, int>, shared_ptr<Arc>> pred;        // key: (vertex_id, time)
  sy_map<tuple<int, int>, double> dist_minus_big_value; // key: (vertex_id, time)

  std::priority_queue<vd, vector<vd>, cmp_vd> queue;
  sy_map<tuple<int, int>, bool> vt_visited; // default is false

  /* Dijkstra search */
  dist_minus_big_value[make_tuple(goal->id, 0)] = -SPP_DIST_BIG_VALUE;
  queue.push(vd(goal, dist_minus_big_value[make_tuple(goal->id, 0)]));

  while (!queue.empty())
  {
    auto top = queue.top();
    queue.pop();
    auto u = top.v;
    if (vt_visited[make_tuple(u->id, 0)])
      continue;
    vt_visited[make_tuple(u->id, 0)] = true;

    plain_distance[make_tuple(u, goal)] = dist_minus_big_value[make_tuple(u->id, 0)] + SPP_DIST_BIG_VALUE;

    if (u == vertex &&
        dist_minus_big_value[make_tuple(u->id, 0)] != 0)
    {
      break;
    }

    if (dist_minus_big_value[make_tuple(u->id, 0)] != 0)
    {
      for (auto &arc : graph.get_outout_arcs(u))
      {
        auto &v = arc->to;
        if (dist_minus_big_value[make_tuple(u->id, 0)] + arc->cost < dist_minus_big_value[make_tuple(v->id, 0)])
        {
          dist_minus_big_value[make_tuple(v->id, 0)] = dist_minus_big_value[make_tuple(u->id, 0)] + arc->cost;
          queue.push(vd(v, dist_minus_big_value[make_tuple(v->id, 0)]));
        }
      }
    }
  }

  return dist_minus_big_value[make_tuple(vertex->id, 0)] + SPP_DIST_BIG_VALUE;
}

double heur(Graph &graph, shared_ptr<Vertex> vertex, shared_ptr<Vertex> goal)
{
  if (vertex->type == Vertex::Type::DUMMY_DESTINATION && vertex->agent_id != goal->agent_id)
    return SPP_DIST_BIG_VALUE;

  if (!plain_distance.contains(make_tuple(vertex, goal)))
  {
    plain_distance[make_tuple(vertex, goal)] = solve_spp_dijkstra(graph, vertex, goal);
  }

  return plain_distance[make_tuple(vertex, goal)];
}

double heur_b(Graph &graph, shared_ptr<Vertex> vertex, shared_ptr<Vertex> goal)
{
  if (vertex->type == Vertex::Type::DUMMY_ORIGIN && vertex->agent_id != goal->agent_id)
    return SPP_DIST_BIG_VALUE;

  if (!plain_distance.contains(make_tuple(vertex, goal)))
  {
    plain_distance[make_tuple(vertex, goal)] = solve_spp_dijkstra_b(graph, vertex, goal);
    if (plain_distance[make_tuple(vertex, goal)] == SPP_DIST_BIG_VALUE)
    println("{} {}  /  {} {}  {}", vertex->name(), vertex->coord.name(), goal->name(), goal->coord.name(), plain_distance[make_tuple(vertex, goal)]);
  }

  return plain_distance[make_tuple(vertex, goal)];
}

tuple<
    vector<shared_ptr<Vertex>>, double, double>
solve_spp_teg_astar(Graph &graph,
                    int agent_id,
                    int num_agent,
                    int ts,
                    shared_ptr<MasterResult> mas_result,
                    ErasedArcs &erased_arcs)
{
  sy_map<tuple<shared_ptr<Vertex>, int>, bool> vt_occupied;
  return solve_spp_teg_astar(graph,
                             agent_id,
                             num_agent,
                             ts,
                             mas_result,
                             erased_arcs,
                             vt_occupied);
}

tuple<
    vector<shared_ptr<Vertex>>, double, double>
solve_spp_teg_astar(Graph &graph,
                    int agent_id,
                    int num_agent,
                    int ts,
                    set<constraint_old_lag> ctrs,
                    shared_ptr<MasterResult> mas_result,
                    sy_map<tuple<shared_ptr<Vertex>, int>, bool> &vt_occupied)
{
  // generate erased_arcs according to ctrs
  ErasedArcs erased_arcs;
  for (auto ctr : ctrs)
  {
    if (ctr.i == agent_id)
    {
      for (auto out_arc : graph.get_out_arcs(ctr.v))
      {
        erased_arcs[make_tuple(out_arc, ctr.t)] = true;
      }
    }
  }

  return solve_spp_teg_astar(graph,
                             agent_id,
                             num_agent,
                             ts,
                             mas_result,
                             erased_arcs,
                             vt_occupied);
}

// Return (path, reduced_cost, cost)
// where path is vector of shared_ptc<Vertex> and contain both dummy origin and destination vertices
tuple<
    vector<shared_ptr<Vertex>>, double, double>
solve_spp_teg_astar(Graph &graph,
                    int agent_id,
                    int num_agent,
                    int ts,
                    shared_ptr<MasterResult> mas_result,
                    ErasedArcs &erased_arcs,
                    sy_map<tuple<shared_ptr<Vertex>, int>, bool> &vt_occupied)
{
  auto o = graph.get_dummy_vertex(Vertex::Type::DUMMY_ORIGIN, agent_id);
  auto d = graph.get_dummy_vertex(Vertex::Type::DUMMY_DESTINATION, agent_id);

  double global_altered_cost = get_global_altered_cost(agent_id, mas_result);

  sy_map<tuple<int, int>, shared_ptr<Arc>> pred;        // key: (vertex_id, time)
  sy_map<tuple<int, int>, double> dist_minus_big_value; // key: (vertex_id, time)

  AstarQueue vertices_to_check;

  /* A* search */
  int t = -1;
  pred[make_tuple(o->id, t)] = nullptr;
  dist_minus_big_value[make_tuple(o->id, t)] = -SPP_DIST_BIG_VALUE;
  double g = dist_minus_big_value[make_tuple(o->id, t)];
  double f = g + heur(graph, o, d);
  vertices_to_check.push(AstarLabel(o, t, g, f));

  while (!vertices_to_check.empty())
  {
    auto top = vertices_to_check.top();
    if (vertices_to_check.empty())
      break;
    vertices_to_check.pop();
    auto &u = top.vertex;
    t = top.time;

    if (u == d)
      break;

    if (dist_minus_big_value[make_tuple(u->id, t)] != 0)
    {
      for (auto &arc : graph.get_out_arcs(u))
      {
        auto &v = arc->to;
        if (erased_arcs.contains(make_tuple(arc, t)) || vt_occupied.contains(make_tuple(v, t + 1)) || vertices_to_check.is_poped(v, t + 1))
          continue;

        double tentative_new_g = dist_minus_big_value[make_tuple(u->id, t)] + arc->cost + get_local_altered_cost(agent_id, mas_result, arc, t);
        if (tentative_new_g < dist_minus_big_value[make_tuple(v->id, t + 1)])
        {
          pred[make_tuple(v->id, t + 1)] = arc;
          dist_minus_big_value[make_tuple(v->id, t + 1)] = tentative_new_g;
          if (mas_result == nullptr)
          {
            vertices_to_check.push(
                AstarLabel(v,
                           t + 1,
                           tentative_new_g,
                           tentative_new_g + heur(graph, v, d),
                           false));
          }
          else
          {
            vertices_to_check.push(
                AstarLabel(v,
                           t + 1,
                           tentative_new_g,
                           tentative_new_g + heur(graph, v, d),
                           mas_result->reserved.contains(make_tuple(v, t + 1))));
          }
        }
      }
    }
  }

  /* Reconstruct path */
  vector<shared_ptr<Vertex>> path;
  double reduced_cost = dist_minus_big_value[make_tuple(d->id, t)] + global_altered_cost + SPP_DIST_BIG_VALUE;
  double cost = 0;

  bool path_reconstructed = true;

  path.insert(path.begin(), d);
  auto &a = pred[make_tuple(d->id, t)];
  for (int tt = t; tt >= 0; tt--)
  {
    if (a == nullptr)
    {
      path_reconstructed = false;
      break;
    }
    path.insert(path.begin(), a->from);
    cost += a->cost;
    a = pred[make_tuple(a->from->id, tt - 1)];
  }

  if (!path_reconstructed || t == -1)
    return make_tuple(vector<shared_ptr<Vertex>>(), DUMMY_COLUMN_COST, DUMMY_COLUMN_COST);

  else
    return make_tuple(path, reduced_cost, cost);
}

void get_forward_dist(Graph &graph,
                      int agent_id,
                      shared_ptr<MasterResult> mas_result,
                      ErasedArcs &erased_arcs,
                      sy_map<tuple<int, int>, double> &dist_minus_big_value_real,
                      sy_map<tuple<int, int>, bool> &vertex_is_in_sub,
                      shared_ptr<Vertex> target_v,
                      int target_t)
{
  auto o = graph.get_dummy_vertex(Vertex::Type::DUMMY_ORIGIN, agent_id);

  AstarQueue vertices_to_check;
  sy_map<tuple<int, int>, double> dist_minus_big_value;

  int t = -1;
  dist_minus_big_value[make_tuple(o->id, t)] = -SPP_DIST_BIG_VALUE;
  double g = dist_minus_big_value[make_tuple(o->id, t)];
  vertices_to_check.push(AstarLabel(o, t, g, g));

  while (!vertices_to_check.empty())
  {
    auto top = vertices_to_check.top();
    if (vertices_to_check.empty())
      break;
    vertices_to_check.pop();
    auto &u = top.vertex;
    t = top.time;

    if (vertex_is_in_sub.contains(make_tuple(u->id, t)))
      dist_minus_big_value_real[make_tuple(u->id, t)] = dist_minus_big_value[make_tuple(u->id, t)];

    if (u == target_v && t == target_t)
      break;

    if (t >= target_t)
      continue;

    if (dist_minus_big_value[make_tuple(u->id, t)] != 0)
    {
      for (auto &arc : graph.get_out_arcs(u))
      {
        auto &v = arc->to;
        if (erased_arcs.contains(make_tuple(arc, t)) || vertices_to_check.is_poped(v, t + 1))
          continue;

        double tentative_new_g = dist_minus_big_value[make_tuple(u->id, t)] + arc->cost + get_local_altered_cost(agent_id, mas_result, arc, t);
        if (tentative_new_g < dist_minus_big_value[make_tuple(v->id, t + 1)])
        {
          dist_minus_big_value[make_tuple(v->id, t + 1)] = tentative_new_g;
          vertices_to_check.push(
              AstarLabel(
                  v,
                  t + 1,
                  tentative_new_g,
                  tentative_new_g + heur(graph, v, target_v),
                  !vertex_is_in_sub.contains(make_tuple(v->id, t + 1)),
                  dist_minus_big_value_real.contains(make_tuple(v->id, t + 1))));
        }
      }
    }
  }
}


void get_backward_dist(Graph &graph,
                      int agent_id,
                      shared_ptr<MasterResult> mas_result,
                      ErasedArcs &erased_arcs,
                      sy_map<tuple<int, int>, double> &dist_minus_big_value_real,
                      sy_map<tuple<int, int>, bool> &vertex_is_in_sub,
                      shared_ptr<Vertex> target_v,
                      int target_t,
                      int tau)
{
  auto d = graph.get_dummy_vertex(Vertex::Type::DUMMY_DESTINATION, agent_id);

  AstarQueue vertices_to_check;
  sy_map<tuple<int, int>, double> dist_minus_big_value;

  int t = tau + b_horizon_margin + 1;
  auto real_d = graph.get_inin_arcs(d)[0]->from;

  for (int tt = t; tt >= target_t + heur(graph, target_v, d); tt--)
  {
    dist_minus_big_value[make_tuple(real_d->id, tt)] = -SPP_DIST_BIG_VALUE;
    vertices_to_check.push(AstarLabel(real_d, tt, dist_minus_big_value[make_tuple(real_d->id, tt)], dist_minus_big_value[make_tuple(real_d->id, tt)]));
  }

  while (!vertices_to_check.empty())
  {
    auto top = vertices_to_check.top();
    if (vertices_to_check.empty())
      break;
    vertices_to_check.pop();
    auto &u = top.vertex;
    t = top.time;

    if (vertex_is_in_sub.contains(make_tuple(u->id, t)))
      dist_minus_big_value_real[make_tuple(u->id, t)] = dist_minus_big_value[make_tuple(u->id, t)];

    if (u == target_v && t == target_t)
      break;

    if (t <= target_t)
      continue;

    if (dist_minus_big_value[make_tuple(u->id, t)] != 0)
    {
      for (auto &arc : graph.get_in_arcs(u))
      {
        auto &v = arc->from;
        if (erased_arcs.contains(make_tuple(arc, t - 1)) || vertices_to_check.is_poped(v, t - 1))
          continue;

        double tentative_new_g = dist_minus_big_value[make_tuple(u->id, t)] + arc->cost + get_local_altered_cost(agent_id, mas_result, arc, t-1);
        if (tentative_new_g < dist_minus_big_value[make_tuple(v->id, t - 1)])
        {
          dist_minus_big_value[make_tuple(v->id, t - 1)] = tentative_new_g;
          vertices_to_check.push(
              AstarLabel(
                  v,
                  t - 1,
                  tentative_new_g,
                  tentative_new_g + heur_b(graph, v, target_v),
                  !vertex_is_in_sub.contains(make_tuple(v->id, t - 1)),
                  dist_minus_big_value_real.contains(make_tuple(v->id, t - 1))));
        }
      }
    }
  }
}

struct vt
{
  shared_ptr<Vertex> v;
  int t;

  vt(shared_ptr<Vertex> v, int t) : v(v), t(t) {}
};

struct cmp_vt
{
  bool operator()(vt l, vt r) { return l.t < r.t; }
};

struct cmp_vt_reverse
{
  bool operator()(vt l, vt r) { return l.t > r.t; }
};
