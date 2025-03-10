#include "cbs_ts.h"
#include "spp_teg.h"
#include <map>
#include <queue>
#include <functional>

#define BIGM 10000000

// MAPFSol::find_first_time_conflict를 참고
collision_t find_collision_t(vector<vector<shared_ptr<Vertex>>> routes, int ts)
{
  int num_agents = routes.size();
  int makespan = 0;
  for (auto &path : routes)
    makespan = makespan < path.size() ? path.size() : makespan;

  sy_map<tuple<shared_ptr<Vertex>, int>, tuple<int, int>> vt_occupied_by; // input: (v,t) / output: (i, t_i)
  for (int t = 0; t <= makespan; t++)
  {
    for (int i = 0; i < num_agents; i++)
    {
      if (t < routes[i].size())
        for (int tt = t; tt <= makespan && tt <= t + ts; tt++)
        {
          if (std::get<0>(vt_occupied_by[make_tuple(routes[i][t], tt)]) == 0)
            vt_occupied_by[make_tuple(routes[i][t], tt)] = make_tuple(i + 1, t);
          else if (std::get<0>(vt_occupied_by[make_tuple(routes[i][t], tt)]) == i + 1)
            vt_occupied_by[make_tuple(routes[i][t], tt)] = make_tuple(i + 1, t);
          else
            return (collision_t){
                std::get<0>(vt_occupied_by[make_tuple(routes[i][t], tt)]) - 1,
                i,
                std::get<1>(vt_occupied_by[make_tuple(routes[i][t], tt)]) - 1, // dummy_origin으로 시작하기 때문에 1 더 빼줌.
                t - 1,
                routes[i][t]};
        }
    }
  }
  return (collision_t){-1, -1, -1, -1, nullptr};
}


void get_ctrs_t(
    cbs_node_t *n,
    int i,
    sy_map<tuple<shared_ptr<Vertex>, int>, bool> &vt_occupied)
{
  cbs_node_t *temp = n;
  vt_occupied = sy_map<tuple<shared_ptr<Vertex>, int>, bool>();
  debugln("\tctrs collecting node {}  agent_id {}", n->idx, i);

  while (temp->depth > 0)
  {
    if (temp->ctr.i == i)
    {
      debugln("\t{}, {}, {}, {}", temp->ctr.i, temp->ctr.t1, temp->ctr.t2, temp->ctr.v->name());
      for (int tt = temp->ctr.t1; tt <= temp->ctr.t2; tt++)
      {
        vt_occupied[make_tuple(temp->ctr.v, tt)] = true;
      }
    }
    temp = temp->parent;
  };
}

vector<vector<shared_ptr<Vertex>>> copy_routes_t(vector<vector<shared_ptr<Vertex>>> routes)
{
  vector<vector<shared_ptr<Vertex>>> temp;
  for (int i = 0; i < routes.size(); i++)
  {
    temp.push_back(vector<shared_ptr<Vertex>>());
    temp.at(i).assign(routes.at(i).begin(), routes.at(i).end());
  }
  return temp;
}

// solve MAPF-TS instance(arc collision) via CBS
void cbs_ts(shared_ptr<Instance> instance)
{
  int &ts = instance->param.time_spacing;
  int &max_time = instance->param.time_limit;

  std::priority_queue<cbs_node_t *, vector<cbs_node_t *>, cbs_node_t_cmp> q;
  int idx = 0, iter = 0, num_spp = 0, num_node = 1;
  cbs_node_t *n;
  int last_value;

  auto start_time = get_now();
  bool time_out = false;

  ErasedArcs erased_arcs;
  sy_map<tuple<shared_ptr<Vertex>, int>, bool> vt_occupied;

  // set root node of CBS tree
  cbs_node_t *root = new cbs_node_t();
  int root_val = 0;

  for (int i = 0; i < instance->param.num_agents; i++)
  {
    auto [route, reduced_cost, cost] = solve_spp_teg_astar(
        instance->graph,
        i,
        instance->param.num_agents,
        ts,
        nullptr,
        erased_arcs,
        vt_occupied);

    root->routes.push_back(route);
    num_spp++;
    root_val += cost;
  }
  root->value = root_val;
  root->parent = NULL;
  root->ctr = (constraint_t){-1, -1, -1, nullptr};
  root->idx = idx++;
  root->depth = 0;
  q.push(root);
  debugln("[DEBUG] Setup: CBS_TS root;\t{}", get_elapsed(start_time));

  // begin CBS iteration
  while (!q.empty())
  {
    if (get_elapsed(start_time) > max_time)
    {
      time_out = true;
      break;
    }
    n = q.top();
    q.pop();

    collision_t coll = find_collision_t(n->routes, ts);
    if (coll.i >= 0)
    {
    }

    if (coll.i == -1)
    { // collision not-detected
      println("optimal\t{}\t{}\t{}\t{}\t{}\t{}", get_elapsed(start_time), n->value, iter, num_spp, num_node, n->depth);
      return;
    }
    else
    { // collision detected
      if (iter % 10000 == 0 || n->value != last_value){
        println("ongoing\t{}\t{}\t{}\t{}\t{}\t{}", get_elapsed(start_time), n->value, iter, num_spp, num_node, n->depth);
        last_value = n->value;
      }
#ifdef DEBUG
      // report detected collision and corresponding routes
      println("coll i{}, j{}, t_i{}, t_j{}, v{}", coll.i, coll.j, coll.t_i, coll.t_j, coll.v->coord.name());

      print("{}: ", coll.i);
      for (auto v : n->routes.at(coll.i))
      {
        print("{} ", v->name());
      }
      println("");

      print("{}: ", coll.j);
      for (auto v : n->routes.at(coll.j))
      {
        print("{} ", v->name());
      }
      println("");
#endif
      cbs_node_t *r = new cbs_node_t();
      cbs_node_t *l = new cbs_node_t();

      r->parent = n;
      r->idx = idx++;
      r->depth = n->depth + 1;
      r->routes = copy_routes_t(n->routes);
      r->ctr = (constraint_t){coll.i, std::min(coll.t_i, coll.t_j), std::min(coll.t_i, coll.t_j) + ts, coll.v};
      get_ctrs_t(r, coll.i, vt_occupied);
      auto [route, reduced_cost, cost] = solve_spp_teg_astar(
          instance->graph,
          coll.i,
          instance->param.num_agents,
          ts,
          nullptr,
          erased_arcs,
          vt_occupied);
      // println("{}, {}", cost, route.size());
      num_spp += 1;
      r->routes.at(coll.i) = route;
      if (r->routes.at(coll.i).size() == 0)
        r->is_leaf = true;
      r->value = n->value - (n->routes.at(coll.i).size() - 3) + cost;

      l->parent = n;
      l->idx = idx++;
      l->depth = n->depth + 1;
      l->routes = copy_routes_t(n->routes);
      l->ctr = (constraint_t){coll.j, std::min(coll.t_i, coll.t_j), std::min(coll.t_i, coll.t_j) + ts, coll.v};
      get_ctrs_t(l, coll.j, vt_occupied);
      std::tie(route, reduced_cost, cost) = solve_spp_teg_astar(
          instance->graph,
          coll.j,
          instance->param.num_agents,
          ts,
          nullptr,
          erased_arcs,
          vt_occupied);
      num_spp += 1;
      l->routes.at(coll.j) = route;
      if (l->routes.at(coll.j).size() == 0)
        l->is_leaf = true;
      l->value = n->value - (n->routes.at(coll.j).size() - 3) + cost;

      if (!r->is_leaf)
      {
        q.push(r);
        num_node += 1;
      }
      if (!l->is_leaf)
      {
        q.push(l);
        num_node += 1;
      }
      iter++;
    }
  } // end of while

  if (time_out)
  {
    println("time_out\t{}", get_elapsed(start_time), iter, num_spp, num_node);
    return;
  }
  else
  {
    println("infeasible", get_elapsed(start_time), iter, num_spp, num_node);
    return;
  }
}
