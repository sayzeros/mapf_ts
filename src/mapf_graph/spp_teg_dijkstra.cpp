#include "spp_teg.h"
#include "sy_log.h"

#include <queue>

#define SPP_DIST_BIG_VALUE 1e9

struct vtd
{
  shared_ptr<Vertex> v;
  int t;
  double dist;

  vtd(shared_ptr<Vertex> v, int t, double dist) : v(v), t(t), dist(dist) {}
};

struct cmp_vtd
{
  bool operator()(vtd l, vtd r) { return l.dist > r.dist; }
};

tuple<
    vector<shared_ptr<Vertex>>, double, double>
solve_spp_teg_dijkstra(Graph &graph,
                       int agent_id,
                       int num_agent,
                       int ts,
                       shared_ptr<MasterResult> mas_result,
                       ErasedArcs &erased_arcs)
{
  sy_map<tuple<shared_ptr<Vertex>, int>, bool> vt_occupied;
  return solve_spp_teg_dijkstra(graph,
                                agent_id,
                                num_agent,
                                ts,
                                mas_result,
                                erased_arcs,
                                vt_occupied);
}

tuple<
    vector<shared_ptr<Vertex>>, double, double>
solve_spp_teg_dijkstra(Graph &graph,
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

  std::priority_queue<vtd, vector<vtd>, cmp_vtd> queue;
  sy_map<tuple<int, int>, bool> vt_visited; // default is false

  /* Dijkstra search */
  int t = -1;
  pred[make_tuple(o->id, t)] = nullptr;
  dist_minus_big_value[make_tuple(o->id, t)] = -SPP_DIST_BIG_VALUE;
  queue.push(vtd(o, t, dist_minus_big_value[make_tuple(o->id, t)]));

  while (!queue.empty())
  {
    auto top = queue.top();
    queue.pop();
    auto u = top.v;
    if (vt_visited[make_tuple(u->id, top.t)])
      continue;
    vt_visited[make_tuple(u->id, top.t)] = true;
    t = top.t;

    if (u == d &&
        dist_minus_big_value[make_tuple(u->id, t)] != 0)
    {
      break;
    }

    if (dist_minus_big_value[make_tuple(u->id, t)] != 0)
    {
      for (auto &arc : graph.get_out_arcs(u))
      {
        if (erased_arcs.contains(make_tuple(arc, t)))
          continue;
        auto &v = arc->to;
        if (vt_occupied.contains(make_tuple(v, t + 1)))
          continue;

        double local_cost_altered = get_local_altered_cost(agent_id, mas_result, arc, t);

        if (dist_minus_big_value[make_tuple(u->id, t)] + arc->cost + local_cost_altered < dist_minus_big_value[make_tuple(v->id, t + 1)])
        {
          pred[make_tuple(v->id, t + 1)] = arc;
          dist_minus_big_value[make_tuple(v->id, t + 1)] = dist_minus_big_value[make_tuple(u->id, t)] + arc->cost + local_cost_altered;
          queue.push(vtd(v, t + 1, dist_minus_big_value[make_tuple(v->id, t + 1)]));
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