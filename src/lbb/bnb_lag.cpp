// #define DEBUG

#ifdef DEBUG
#define D(x) (x)
#define ND(x) \
  do          \
  {           \
  } while (0)
#include <iostream>
#else
#define D(x) \
  do         \
  {          \
  } while (0)
#define ND(x) (x)
#endif

#include "bnb_lag.h"
#include "spp_teg.h"

#include <cmath>
#include <map>
#include <functional>
#include <utility>

using Ctrs_t = set<constraint_old_lag>;

Ctrs_t get_ctrs_t(bnb_node *n)
{
  Ctrs_t ctrs;
  bnb_node *temp = n;

  while (temp->depth > 0)
  {
    if (temp->ctr.i != -1)
    {
      ctrs.insert(temp->ctr);
    }
    temp = temp->parent;
  }
  return ctrs;
}

collision_t find_first_collision(vector<vector<shared_ptr<Vertex>>> routes, int ts)
{
  // check makespan of routes
  int max_len = 0;
  for (int i = 0; i < routes.size(); i++)
  {
    max_len = max_len > routes.at(i).size() ? max_len : routes.at(i).size();
  }

  // detect collision
  for (int t1 = 0; t1 < max_len; t1++)
  {
    for (int i = 0; i < routes.size(); i++)
    {
      for (int j = i + 1; j < routes.size(); j++)
      {
        shared_ptr<Vertex> temp1, temp2;

        if (t1 < routes.at(i).size())
          temp1 = routes.at(i).at(t1);
        else
          continue;

        for (int tt = -ts; tt <= ts; tt++)
        {
          if (t1 + tt >= 0)
          {
            if (t1 + tt < routes.at(j).size())
              temp2 = routes.at(j).at(t1 + tt);
            else
              continue;
            if (temp1 == temp2)
            {
              return (collision_t){i, j, t1 - 1, t1 + tt - 1, temp1}; // -1 since route is start with dummy_origin
            }
          }
        }
      }
    }
  }
  return (collision_t){-1, -1, -1, -1, nullptr};
}

int update_multiplier(shared_ptr<Instance> instance,
                      vector<vector<shared_ptr<Vertex>>> &routes,
                      vt_values &lamb_v,
                      at_values &lamb_a,
                      ivt_values &lamb_o,
                      vt_values &lamb_r,
                      double it_cost,
                      double lb, // best lb of the node
                      double ub,
                      int lag_it)
{
  Graph &graph = instance->graph;
  int &ts = instance->param.time_spacing;
  string ctrs_string = instance->param.formul_string + instance->param.ctrs_string;

  vt_values cons_v_lhs;
  at_values cons_a_lhs;
  ivt_values cons_o_lhs;
  vt_values cons_r_lhs;

  double alpha_v, alpha_a, alpha_o, alpha_r;

  int n_col = 0;
  if (ctrs_string.find('V') != std::string::npos)
    for (auto &[vt, val] : lamb_v)
    {
      cons_v_lhs[vt] = 0; 
    }
  if (ctrs_string.find('A') != std::string::npos)
    for (auto &[at, val] : lamb_a)
    {
      cons_a_lhs[at] = 0;
    }
  if (ctrs_string.find('O') != std::string::npos)
    for (auto &[ivt, val] : lamb_o)
    {
      cons_o_lhs[ivt] = 0;
    }
  if (ctrs_string.find('R') != std::string::npos)
    for (auto &[vt, val] : lamb_r)
    {
      cons_r_lhs[vt] = 0;
    }

  // calculate inequalities' lhs
  int i = -1;
  for (auto &route : routes)
  {
    i++;
    debug_assert(route.size() > 1);

    int t = -2; 
    shared_ptr<Vertex> last_v = nullptr;
    for (auto v : route)
    {
      t++; 

      // Constraint V
      cons_v_lhs[make_tuple(v->id, t)] += 1;
      // Constraint O
      cons_o_lhs[make_tuple(i, v->id, t)] += 1;
      for (int j = 0; j < routes.size(); j++)
      {
        if (j != i)
          for (int delt = -ts; delt <= 0; delt++)
            cons_o_lhs[make_tuple(j, v->id, t + delt)] += 1 / (double)(ts + 1);
      }
      // Constraint A & R
      if (last_v != nullptr)
      {
        // A
        if (last_v->id < v->id)
          cons_a_lhs[make_tuple(graph.get_arc(last_v, v)->id, t - 1)] += 1;
        else if (last_v->id > v->id && graph.get_arc(v, last_v) != nullptr)
          cons_a_lhs[make_tuple(graph.get_arc(v, last_v)->id, t - 1)] += 1;
        // R
        if (last_v == v)
        {
          cons_r_lhs[make_tuple(v->id, t - ts - 1)] += 1;
        }
        else // last_v != v
        {
          for (int delt = -ts - 1; delt <= -1; delt++)
            cons_r_lhs[make_tuple(last_v->id, t + delt)] += 1;
        }
      }
      last_v = v;
    }
  }

  // calculate alpha_v
  alpha_v = 0;
  // update lamb_v
  if (ctrs_string.find('V') != std::string::npos)
    for (const auto &[vt, count] : cons_v_lhs)
    {
      lamb_v[vt] += alpha_v * (count - 1);
      if (lamb_v[vt] <= 0)
        lamb_v.erase(vt);
      if (count > 1)
        n_col += 1;
    }
  // calculate alpha_a
  if (instance->param.stepsize_str.find("fixed1") != std::string::npos)
    alpha_a = 0.001;
  else if (instance->param.stepsize_str.find("fixed2") != std::string::npos)
    alpha_a = 0.005;
  else if (instance->param.stepsize_str.find("fixed3") != std::string::npos)
    alpha_a = 0.01;
  else if (instance->param.stepsize_str.find("fixed4") != std::string::npos)
    alpha_a = 0.05;
  else if (instance->param.stepsize_str.find("fixed5") != std::string::npos)
    alpha_a = 0.1;
  else if (instance->param.stepsize_str.find("dimini") != std::string::npos)
    alpha_a = 1 / ((double)(lag_it + 1));
  else // if (instance->param.stepsize_str.find("polyak") != std::string::npos)
  {
    if (lag_it < 50)
      alpha_a = 1 / ((double)(lag_it + 1));
    else
    {
      double sum_square_subgrad = 0;
      for (const auto &[vt, count] : cons_a_lhs)
        sum_square_subgrad += (count - 1) * (count - 1);
      double apx_target = lb + (ub - lb) * (1 / std::sqrt((double)(lag_it + 1)));
      if (sum_square_subgrad == 0)
        alpha_a = 0;
      else
        alpha_a = std::min(apx_target - it_cost, 5.0) / sum_square_subgrad;
    }
  }
  // update lamb_a
  if (ctrs_string.find('A') != std::string::npos)
    for (const auto &[at, count] : cons_a_lhs)
    {
      lamb_a[at] += alpha_a * (count - 1);
      if (lamb_a[at] <= 0)
        lamb_a.erase(at);
      if (count > 1)
      {
        n_col += 1;
      }
    }
  // calculate alpha_o
  if (instance->param.stepsize_str.find("fixed1") != std::string::npos)
    alpha_o = 0.001;
  else if (instance->param.stepsize_str.find("fixed2") != std::string::npos)
    alpha_o = 0.005;
  else if (instance->param.stepsize_str.find("fixed3") != std::string::npos)
    alpha_o = 0.01;
  else if (instance->param.stepsize_str.find("fixed4") != std::string::npos)
    alpha_o = 0.05;
  else if (instance->param.stepsize_str.find("fixed5") != std::string::npos)
    alpha_o = 0.1;
  else if (instance->param.stepsize_str.find("dimini") != std::string::npos)
    alpha_o = 1 / ((double)(lag_it + 1));
  else // if (instance->param.stepsize_str.find("polyak") != std::string::npos)
  {
    if (lag_it < 50)
      alpha_o = 1 / ((double)(lag_it + 1));
    else
    {
      double sum_square_subgrad = 0;
      for (const auto &[vt, count] : cons_o_lhs)
        sum_square_subgrad += (count - 1) * (count - 1);
      double apx_target = lb + (ub - lb) * (1 / std::sqrt((double)(lag_it + 1)));
      if (sum_square_subgrad == 0)
        alpha_o = 0;
      else
        alpha_o = std::min(apx_target - it_cost, 5.0) / sum_square_subgrad;
      if (alpha_o < 0)
      {
        println("o {} {} {}", alpha_a, apx_target, it_cost);
        exit(0);
      }
    }
  }
  // update lamb_o
  if (ctrs_string.find('O') != std::string::npos)
    for (const auto &[ivt, count] : cons_o_lhs)
    {
      lamb_o[ivt] += alpha_o * (count - 1);
      if (lamb_o[ivt] <= 0)
        lamb_o.erase(ivt);
      if (count > 1)
      {
        n_col += 1;
      }
    }
  // calculate alph_r
  if (instance->param.stepsize_str.find("fixed1") != std::string::npos)
    alpha_r = 0.001;
  else if (instance->param.stepsize_str.find("fixed2") != std::string::npos)
    alpha_r = 0.005;
  else if (instance->param.stepsize_str.find("fixed3") != std::string::npos)
    alpha_r = 0.01;
  else if (instance->param.stepsize_str.find("fixed4") != std::string::npos)
    alpha_r = 0.05;
  else if (instance->param.stepsize_str.find("fixed5") != std::string::npos)
    alpha_r = 0.1;
  else if (instance->param.stepsize_str.find("dimini") != std::string::npos)
    alpha_r = 1 / ((double)(lag_it + 1));
  else // if (instance->param.stepsize_str.find("polyak") != std::string::npos)
  {
    if (lag_it < 50)
      alpha_r = 1 / ((double)(lag_it + 1));
    else
    {
      double sum_square_subgrad = 0;
      for (const auto &[vt, count] : cons_r_lhs)
        sum_square_subgrad += (count - 1) * (count - 1);
      double apx_target = lb + (ub - lb) * (1 / std::sqrt((double)(lag_it + 1)));
      if (sum_square_subgrad == 0)
        alpha_r = 0;
      else
        alpha_r = std::min(apx_target - it_cost, 5.0) / sum_square_subgrad;
    }
  }
  // update lamb_r TODO:
  if (ctrs_string.find('R') != std::string::npos)
    for (const auto &[vt, count] : cons_r_lhs)
    {
      lamb_r[vt] += alpha_r * (count - 1);
      if (lamb_r[vt] <= 0)
        lamb_r.erase(vt);
      if (count > 1)
      {
        n_col += 1;
      }
    }

  return n_col;
}

tuple<vector<vector<shared_ptr<Vertex>>>, double> fixing_heuristic(vector<vector<shared_ptr<Vertex>>> &routes, shared_ptr<Instance> instance, bnb_node *node)
{
  // following thesis using option 8 (not using original path, original path length decreasing order, and using updated cost by LR iter)
  vector<vector<shared_ptr<Vertex>>> temp_routes;
  temp_routes.resize(routes.size());

  // generate re-routing order in decreasing length of the original routes
  vector<int> order;
  int idx = -1;
  for (auto &route : routes)
  {
    idx++;
    if (idx == 0)
    {
      order.push_back(idx);
      continue;
    }
    for (int a = 0; a < idx; a++)
    {
      if (route.size() < routes[order[a]].size())
      {
        order.insert(order.begin() + a, idx);
        break;
      }
    }
    if (order.size() < idx + 1)
      order.push_back(idx);
  }

  // generate new routes
  double total_cost = 0;
  sy_map<tuple<shared_ptr<Vertex>, int>, bool> reserved;
  auto mas_result = make_shared<MasterResult>();
  // reflect lamb_v into mas_result
  for (auto &[key, value] : node->lamb_v)
  {
    auto &[v, t] = key;
    for (auto &a : instance->graph.get_in_arcs(instance->graph.get_vertex(v)))
    {
      mas_result->at_cost_altered[make_tuple(a, t - 1)] += value;
    }
  }
  // lamb_a into mas_result
  for (auto &[key, value] : node->lamb_a)
  {
    auto &[a_id, t] = key;
    auto a = instance->graph.get_arc(a_id);
    mas_result->at_cost_altered[make_tuple(a, t)] += value;
    mas_result->at_cost_altered[make_tuple(instance->graph.get_arc(a->to, a->from), t)] += value;
  }
  // reflect lamb_o into mas_result
  for (auto &[key, val] : node->lamb_o)
  {
    auto &[i, v_id, t] = key;
    auto v = instance->graph.get_vertex(v_id);
    for (auto &a : instance->graph.get_in_arcs(v))
    {
      mas_result->i_at_cost_altered[i][make_tuple(a, t - 1)] += val;
      for (int j = 0; j < instance->param.num_agents; j++)
        if (j != i)
          for (int delt = -instance->param.time_spacing; delt <= 0; delt++)
            mas_result->i_at_cost_altered[j][make_tuple(a, t - 1 - delt)] += val / (double)(instance->param.time_spacing + 1);
    }
  }
  // reflect lamb_r into mas_result 
  for (auto &[key, val] : node->lamb_r)
  {
    auto &[v_id, t] = key;
    auto v = instance->graph.get_vertex(v_id);
    mas_result->at_cost_altered[make_tuple(instance->graph.get_arc(v, v), t + instance->param.time_spacing)] += val;
    for (auto &a : instance->graph.get_outout_arcs(v))
      for (int delt = 0; delt <= instance->param.time_spacing; delt++)
        mas_result->at_cost_altered[make_tuple(a, t + delt)] += val;
  }

  for (int i : order)
  {
    double cost, reduced_cost;
    tie(temp_routes[i], reduced_cost, cost) = solve_spp_teg_astar(
        instance->graph,
        i,
        instance->param.num_agents,
        instance->param.time_spacing,
        get_ctrs_t(node),
        mas_result,
        reserved);

    total_cost += cost;

    int t = -2;
    for (auto &v : temp_routes[i])
    {
      t++;
      for (int delt = -instance->param.time_spacing; delt <= instance->param.time_spacing; delt++)
        reserved[make_tuple(v, t + delt)] = true;
    }
  }

  return tie(temp_routes, total_cost);
}

string before_return_update_coll(BnB_Tree &tree, bnb_node *node, vector<vector<shared_ptr<Vertex>>> &routes, shared_ptr<Instance> instance, int &num_spp, string trt)
{
  node->coll = find_first_collision(routes, instance->param.time_spacing);
  if (node->coll.i == -1)
  {
    sy_map<tuple<shared_ptr<Vertex>, int>, bool> not_reserved;
    shared_ptr<MasterResult> not_mas_result;
    // find collision for branching
    routes.clear();
    for (int i = 0; i < instance->param.num_agents; i++)
    {
      vector<shared_ptr<Vertex>> route;
      double route_cost, reduced_cost;

      tie(route, reduced_cost, route_cost) = solve_spp_teg_astar(
          instance->graph,
          i,
          instance->param.num_agents,
          instance->param.time_spacing,
          get_ctrs_t(node),
          not_mas_result,
          not_reserved);

      routes.push_back(route);
      num_spp++;
    }
    node->coll = find_first_collision(routes, instance->param.time_spacing);
    if (node->coll.i == -1)
    {
      tree.ub_routes = routes;
      return "optimal";
    }
  }
  return trt;
}

// Lagrangian for root node
string solve_LR(BnB_Tree &tree, bnb_node *node, int max_lag_iter, int &num_spp, shared_ptr<Instance> instance, time_point start_time)
{
  string ctrs_string = instance->param.formul_string + instance->param.ctrs_string;

  double last_lb = node->lb;
  double eps_change = 0.01;
  int changed_since = 0;
  int not_changing_criteria = 100; 

  int n_col, lag_it;
  string status = "proceed";

  // Lagrangian iteration starts!
  vector<vector<shared_ptr<Vertex>>> routes;
  for (lag_it = 0; lag_it <= max_lag_iter; lag_it++)
  {
    // iteration cost
    double it_cost = 0;
    if (ctrs_string.find('V') != std::string::npos)
      for (const auto &[vt, val] : node->lamb_v)
        it_cost -= val;
    if (ctrs_string.find('A') != std::string::npos)
      for (const auto &[vt, val] : node->lamb_a)
        it_cost -= val;
    if (ctrs_string.find('O') != std::string::npos)
      for (const auto &[vt, val] : node->lamb_o)
        it_cost -= val;
    if (ctrs_string.find('R') != std::string::npos)
      for (const auto &[vt, val] : node->lamb_r)
        it_cost -= val;

    //// re-calculate routes for agents
    sy_map<tuple<shared_ptr<Vertex>, int>, bool> not_reserved;
    auto mas_result = make_shared<MasterResult>();
    // reflect lamb_v into mas_result
    for (auto &[key, value] : node->lamb_v)
    {
      auto &[v, t] = key;
      for (auto a : instance->graph.get_in_arcs(instance->graph.get_vertex(v)))
      {
        mas_result->at_cost_altered[make_tuple(a, t - 1)] += value;
      }
    }
    // reflect lamb_a into mas_result
    for (auto &[key, value] : node->lamb_a)
    {
      auto &[a_id, t] = key;
      auto a = instance->graph.get_arc(a_id);
      mas_result->at_cost_altered[make_tuple(a, t)] += value;
      mas_result->at_cost_altered[make_tuple(instance->graph.get_arc(a->to, a->from), t)] += value;
    }
    // reflect lamb_o into mas_result
    for (auto &[key, val] : node->lamb_o)
    {
      auto &[i, v_id, t] = key;
      auto v = instance->graph.get_vertex(v_id);
      for (auto &a : instance->graph.get_in_arcs(v))
      {
        mas_result->i_at_cost_altered[i][make_tuple(a, t - 1)] += val;
        for (int j = 0; j < instance->param.num_agents; j++)
          if (j != i)
            for (int delt = -instance->param.time_spacing; delt <= 0; delt++)
              mas_result->i_at_cost_altered[j][make_tuple(a, t - 1 - delt)] += val / (double)(instance->param.time_spacing + 1);
      }
    }
    // reflect lamb_r into mas_result
    for (auto &[key, val] : node->lamb_r)
    {
      auto &[v_id, t] = key;
      auto v = instance->graph.get_vertex(v_id);
      mas_result->at_cost_altered[make_tuple(instance->graph.get_arc(v, v), t + instance->param.time_spacing)] += val;
      for (auto &a : instance->graph.get_outout_arcs(v))
        for (int delt = 0; delt <= instance->param.time_spacing; delt++)
          mas_result->at_cost_altered[make_tuple(a, t + delt)] += val;
    }

    routes.clear();
    for (int i = 0; i < instance->param.num_agents; i++)
    {
      vector<shared_ptr<Vertex>> route;
      double route_cost, reduced_cost;

      tie(route, reduced_cost, route_cost) = solve_spp_teg_astar(
          instance->graph,
          i,
          instance->param.num_agents,
          instance->param.time_spacing,
          get_ctrs_t(node),
          mas_result,
          not_reserved);
          
      routes.push_back(route);
      num_spp++;
      it_cost += reduced_cost;
    }
    
    //// update the lower bound
    if (node->lb < it_cost)
    {
      node->lb = it_cost;
    }
    if (node->lb > tree.ub || it_cost > 1e6)
    {
      return "infeasible";
    }

    if (abs(last_lb - node->lb) < eps_change)
    {
      changed_since += 1;
    }
    else
    {
      changed_since = 0;
      last_lb = node->lb;
    }

    //// find violated constraints and calculate Lagrangian multiplier
    n_col = update_multiplier(instance, routes, node->lamb_v, node->lamb_a, node->lamb_o, node->lamb_r, it_cost, node->lb, tree.ub, lag_it);

    if (n_col == 0)
    { // somehow current routes are feasible!
      double ub = 0;
      for (auto &route : routes)
      {
        ub += route.size() - 1;
      }
      if (tree.ub > ub)
      {
        tree.ub = ub;
        tree.ub_routes = routes;
      }
    }
    else if (changed_since == 0)
    {
      // only apply heuristic when lower bound change more than eps
      double heur_cost;
      vector<vector<shared_ptr<Vertex>>> heur_routes;
      tie(heur_routes, heur_cost) = fixing_heuristic(routes, instance, node);
      // update upperbound if it would
      if (tree.ub > heur_cost)
      {
        tree.ub = heur_cost;
        tree.ub_routes = heur_routes;
      }
    }

    if (changed_since == 0 || lag_it % 1000 == 0)
      println("    lag_it {:3}:  time {:.4f}  lb {:.6f}  ub {}  #col {}  not_changed_for {}", lag_it, get_elapsed(start_time), node->lb, tree.ub, n_col, changed_since);

    if (changed_since >= not_changing_criteria)
    {
      // only apply heuristic when lower bound change more than eps
      double heur_cost;
      vector<vector<shared_ptr<Vertex>>> heur_routes;
      tie(heur_routes, heur_cost) = fixing_heuristic(routes, instance, node);
      // update upperbound if it would
      if (tree.ub > heur_cost)
      {
        tree.ub = heur_cost;
        tree.ub_routes = heur_routes;
      }
      status = "proceed";
      break;
    }
    if (get_elapsed(start_time) > instance->param.time_limit)
    {
      status = "time-out";
      break;
    }
    if (tree.ub - ceil(node->lb) < 1)
    {
      status = "optimal";
      break;
    }
  }
  println("    lag_it {:3}:  time {:.4f}  lb {:.6f}  ub {}  #col {}  not_changed_for {}", lag_it, get_elapsed(start_time), node->lb, tree.ub, n_col, changed_since);
  return before_return_update_coll(tree, node, routes, instance, num_spp, status);
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////   Main Part of This Source Code   //////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
void bnb_lag(shared_ptr<Instance> instance)
{
  // basic settings
  auto start_time = get_now();
#ifdef DEBUG
  println("Starting [{}]", get_elapsed(start_time));
#endif

  // algorithm setting
  int max_lag_iter = 5000; // FIXME:

  // algorithm trackers
  int iter = 0;
  int num_spp = 0;
  string status;
  double global_lb;

  // set bnb_tree and root node
  BnB_Tree tree;
  bnb_node *root = tree.generate_root();

  status = solve_LR(tree, root, max_lag_iter, num_spp, instance, start_time);

  // return; // FIXME: node 1개만 풀게 함.
  tree.push(root);

  if (status != "proceed")
  {
    tree.check_out(root);
    global_lb = root->lb;
    tree.pop();
  }

  //////////////////////////
  // begin BB iteration ///
  ////////////////////////
  while (!tree.empty())
  {
    if (get_elapsed(start_time) > instance->param.time_limit)
    {
      status = "time_out";
      break;
    }
    bnb_node *n = tree.top();
    global_lb = n->lb;
    tree.pop();
    
    if (tree.ub - global_lb < 1)
      break;

    // #ifdef DEBUG
    if (n->parent != NULL)
      println("bnb_it {:2}: time {:.6f}  lb {:.6f}  ub {}  depth {}  node {:2}  child_of {:2}", iter, get_elapsed(start_time), n->lb, tree.ub, n->depth, n->idx, n->parent->idx);
    else
      println("bnb_it {:2}: time {:.6f}  lb {:.6f}  ub {}  depth {}  node {:2}", iter, get_elapsed(start_time), n->lb, tree.ub, n->depth, n->idx);
    // #endif
    release_assert(get_ctrs_t(n).size() == n->depth);


    bnb_node *l;
    bnb_node *r;
    std::tie(l, r) = tree.get_branches(n);

    // solve left child
    if (l->ctr.v == instance->graph.get_dummy_vertex(Vertex::Type::DUMMY_ORIGIN, l->ctr.i) && (l->ctr.t <= instance->param.time_spacing || l->ctr.v == instance->graph.get_dummy_vertex(Vertex::Type::DUMMY_DESTINATION, l->ctr.i)))
    {
      l->is_leaf = true;
      tree.check_out(l);
    }
    else
    {
      // #ifdef DEBUG
      println("  l {} of node {}", l->idx, n->idx);
      // #endif
      status = solve_LR(tree, l, max_lag_iter, num_spp, instance, start_time);
      if (status == "proceed")
      {
        tree.push(l);
      }
      // #ifdef DEBUG
      println("  l {}\n", status);
      // #endif
    }

    // solve right child
    if (r->ctr.v == instance->graph.get_dummy_vertex(Vertex::Type::DUMMY_ORIGIN, r->ctr.i) && (r->ctr.t <= instance->param.time_spacing || r->ctr.v == instance->graph.get_dummy_vertex(Vertex::Type::DUMMY_DESTINATION, r->ctr.i)))
    {
      r->is_leaf = true;
      tree.check_out(r);
    }
    else
    {
      // #ifdef DEBUG
      println("  r {} of node {}", r->idx, n->idx);
      // #endif
      status = solve_LR(tree, r, max_lag_iter, num_spp, instance, start_time);
      if (status == "proceed")
      {
        tree.push(r);
      }
      // #ifdef DEBUG
      println("  r {}\n", status);
      // #endif
    }

    iter++;
    tree.check_out(n);

    if (status == "time-out")
    {
      break;
    }
  }
  if (tree.ub - global_lb < 1)
    status = "optimal";

  println("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", "Time", "Status", "lb", "#Node", "BestBound", "RootBound", "#Visited", "#Spp");
  println("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", get_elapsed(start_time), status, global_lb, tree.num_node, tree.ub, root->lb, tree.num_visited, num_spp);
}