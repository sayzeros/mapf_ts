#include "bpc_tree.h"

/* Constructor which will be called to generate & solve root node */
BPCTree::BPCTree(shared_ptr<Instance> instance) : instance(instance)
{
  tree_start_time = get_now();
  bound_updated = false;
  status = Status::ONGOING;

  // Column pool
  pool = make_shared<ColumnPool>();

  // Time records
  elapsed_time = 0.0;
  total_time_on_mp = 0.0;
  total_time_on_sp = 0.0;
  total_time_on_sepa = 0.0;
  total_time_on_branching = 0.0;
  total_time_on_heur = 0.0;
  total_time_on_node = 0.0;

  avg_time_on_erasing = 0.0;
  avg_time_on_sp = 0.0;

  num_spp = 0;

  num_vc_ctrs = 0, 
  num_ac_ctrs = 0;
  num_tc_ctrs = 0;
  num_wc_ctrs = 0;
  num_rc_ctrs = 0;
  num_oc_ctrs = 0;
  num_cc_ctrs = 0;
  num_1c_ctrs = 0;
  num_2c_ctrs = 0;
  num_3c_ctrs = 0;
  num_4c_ctrs = 0;

  // Priority queue that stores unexplored nodes (priority is based on the parent's lower bound)
  unexplored_nodes = NodeQueue();

  // Initialize statistics
  max_depth = 0;
  num_nodes_generated = 0;
  ub = DOUBLE_NULL;
  lb = DOUBLE_NULL;
  node_attaining_ub = nullptr;
  gap_at_root = DOUBLE_NULL;
  gap = DOUBLE_NULL;

  // Initialize solvers
  mas_solver = make_shared<MasterSolver>(instance);
  sub_solver = make_shared<SubSolver>(instance);

  // Generate root
  auto root = make_shared<BPCNode>(instance, pool, mas_solver, sub_solver);
  unexplored_nodes.push(root);
  num_nodes_generated++;

#ifdef SHOW_BBNODE_RESULT
  elapsed_time = get_elapsed(tree_start_time);
  println("Tree: Initial-set-up-done  time {:.3f}", elapsed_time); // 0 initial spaces
#endif
}

void BPCTree::solve_root(bool called_by_user)
{
  auto root_start_time = get_now();
  auto root = unexplored_nodes.top();
  root->solve();
  elapsed_time = get_elapsed(root_start_time);
  update_stats(root);
  status = root->is_feasible() ? Status::LP_OPTIMAL : Status::INFEASIBLE;
  if (elapsed_time >= instance->param.time_limit)
    status = Status::TIME_LIMIT_REACHED;

  if (called_by_user)
  {
    print_result();
  }
}

void BPCTree::explore_tree()
{
  shared_ptr<BPCNode> current_node;
  while (!unexplored_nodes.empty())
  {
    current_node = unexplored_nodes.top();
    unexplored_nodes.pop();

    if (lb != current_node->parent_lb)
    {
      lb = current_node->parent_lb;
      bound_updated = true;
    }
    if (current_node->id > 1 && ub <= std::ceil(lb - INTEGRALITY_EPS))
      break; // current ub is optimal

    auto node_solve_start_time = get_now();
    current_node->solve();
    total_time_on_node += get_elapsed(node_solve_start_time);

    bool node_feasible = current_node->is_feasible();
    if (!node_feasible)
    { // prune by infeasibility
      update_stats(current_node);
      continue;
    }
    bool node_fractional = current_node->has_fractional_solution();
    if (!node_fractional)
    {
      if (current_node->lb <= ub)
      {
        ub = current_node->lb;
        bound_updated = true;
      }
    }
    else if (node_fractional)
    {
      primal_heuristic(current_node);
    }
    update_stats(current_node);
    if (ub <= std::ceil(lb - INTEGRALITY_EPS))
      break; // current ub is optimal

    bool node_suboptimal = std::ceil(current_node->lb - INTEGRALITY_EPS) >= ub;
#ifdef SHOW_BBNODE_RESULT
#ifndef SHOW_BBNODE_RESULT_MINIMAL
    current_node->print_result(unexplored_nodes.size(), node_feasible, node_suboptimal, node_fractional, ub);
#endif
#endif
    if (node_suboptimal)
    { // prune by sub-optimality
      continue;
    }

    if (node_fractional)
    {
      auto branching_start_time = get_now();
      branch(current_node); // while branching, give child node appropriate id and increase num_nodes_generated
      total_time_on_branching += get_elapsed(branching_start_time);
    }

#ifdef SHOW_BBNODE_RESULT_MINIMAL
    if (current_node->id == 1)
    {
      println("  Node-{}: time {:.4f}  lb -  ub -  %gap -", BPCNode::num_explored, get_elapsed(tree_start_time));
      last_lb = current_node->parent_lb;
      last_ub = ub;
    }
    else if (last_lb != current_node->parent_lb || last_ub != ub)
    {
      println("  Node-{}: time {:.4f}  lb {:.5f}  ub {}  %gap {:.3f}", BPCNode::num_explored, get_elapsed(tree_start_time), current_node->parent_lb, (int)ub, (ub - current_node->parent_lb) / ub * 100);
      last_lb = current_node->parent_lb;
      last_ub = ub;
    }
#endif

    elapsed_time = get_elapsed(tree_start_time);

    if (elapsed_time >= instance->param.time_limit)
    {
      println("Tree: time-limit-reached {}sec", instance->param.time_limit); // 0 initial spaces
      status = Status::TIME_LIMIT_REACHED;
      break;
    }
  }
  /* Check and Update status */
  if (status != Status::TIME_LIMIT_REACHED)
  {
    if (ub != DOUBLE_NULL)
    {
      status = Status::OPTIMAL_FOUNDED;
      lb = ub;
      gap = 0;
    }
    else
    {
      status = Status::INFEASIBLE;
    }
  }
  string status_string = get_status_string();

#ifdef SHOW_BBNODE_RESULT_MINIMAL
  println("  Node-{}: time {:.4f}  lb {:.5f}  ub {}  %gap {:.3f}", BPCNode::num_explored, get_elapsed(tree_start_time), current_node->parent_lb, (int)ub, (ub - current_node->parent_lb) / ub * 100);
#endif
  elapsed_time = get_elapsed(tree_start_time);

  print_result();
}

void BPCTree::primal_heuristic(shared_ptr<BPCNode> current_node)
{
  auto heur_start = get_now();

  auto mapf_sol = current_node->primal_heuristic();
  if (mapf_sol.obj_val < ub)
  {
    ub = mapf_sol.obj_val;
    bound_updated = true;
  }

  total_time_on_heur += get_elapsed(heur_start);
}

void BPCTree::branch(shared_ptr<BPCNode> node)
{
  branch_on_forbid_arcs(node); 
}

class ArcTimeCompare 
{
public:
  bool operator()(const tuple<shared_ptr<Arc>, int> &at1, const tuple<shared_ptr<Arc>, int> &at2) const
  {
    return (std::get<1>(at1) > std::get<1>(at2));
  }
};

using ArcTimeQueue = std::priority_queue<
    tuple<shared_ptr<Arc>, int>,
    vector<tuple<shared_ptr<Arc>, int>>,
    ArcTimeCompare>;

bool BPCTree::branch_on_forbid_arcs(shared_ptr<BPCNode> node)
{
  static int branched_till = 0;
  sy_map<int, ArcTimeQueue> arcs_topo_sorted;                         // by agent
  sy_map<int, sy_map<tuple<shared_ptr<Arc>, int>, double>> arc_value; // by agent

  for (auto &[col, val] : node->base_columns)
  {
    int time = -2;
    for (auto &v : col->path)
    {
      time += 1;
      if (time == col->path.size() - 2)
        break;
      auto arc = instance->graph.get_arc(v, col->path[time + 2]);
      if (arc_value[col->agent_id][make_tuple(arc, time)] == 0)
      {
        arcs_topo_sorted[col->agent_id].push(make_tuple(arc, time));
      }
      arc_value[col->agent_id][make_tuple(arc, time)] += val;
    }
  }

  // Detect first diverging vertex
  int fractional = -1;
  tuple<shared_ptr<Vertex>, int> diverging_vt;
  for (int agent_id = 0; agent_id < instance->param.num_agents; agent_id++)
  {
    if (agent_id - 1 < branched_till)
      continue;
    while (!arcs_topo_sorted[agent_id].empty())
    {
      auto at = arcs_topo_sorted[agent_id].top();
      arcs_topo_sorted[agent_id].pop();
      auto val = arc_value[agent_id][at];

      auto [a, t] = at;
      if (val < 1 - INTEGRALITY_EPS)
      {
        fractional = agent_id;
        diverging_vt = make_tuple(std::get<0>(at)->from, std::get<1>(at));
        break;
      }
    }
    if (fractional != -1)
      break;
  }
  if (fractional == -1)
  {
    for (int agent_id = 0; agent_id < instance->param.num_agents; agent_id++)
    {
      if (agent_id > branched_till)
        break;
      while (!arcs_topo_sorted[agent_id].empty())
      {
        auto at = arcs_topo_sorted[agent_id].top();
        arcs_topo_sorted[agent_id].pop();
        auto val = arc_value[agent_id][at];
        if (val < 1 - INTEGRALITY_EPS)
        {
          fractional = agent_id;
          diverging_vt = make_tuple(std::get<0>(at)->from, std::get<1>(at));
          break;
        }
      }
      if (fractional != -1)
        break;
    }
  }
  branched_till += 1;
  if (branched_till >= instance->param.num_agents)
    branched_till -= instance->param.num_agents;

  // Divide outgoing arcs into two set
  sy_map<tuple<shared_ptr<Arc>, int>, bool> forbided_1, forbided_2;
  double cummul_val = 0;
  for (auto &arc : instance->graph.get_out_arcs(std::get<0>(diverging_vt)))
  {
    auto &t = std::get<1>(diverging_vt);
    if ((cummul_val == 0) ||
        (!(arc_value[fractional][make_tuple(arc, t)] > 0 && forbided_2.size() == 0) &&
         (cummul_val + arc_value[fractional][make_tuple(arc, t)] < 0.5 - INTEGRALITY_EPS)))
    {
      forbided_1[make_tuple(arc, t)] = true;
      cummul_val += arc_value[fractional][make_tuple(arc, t)];
    }
    else
    {
      forbided_2[make_tuple(arc, t)] = true;
    }
  }

  shared_ptr<BranchingRule> forbid_arcs_1 = make_shared<ForbidArcs>(fractional, forbided_1);
  shared_ptr<BranchingRule> forbid_arcs_2 = make_shared<ForbidArcs>(fractional, forbided_2);

  unexplored_nodes.push(make_shared<BPCNode>(*node, forbid_arcs_1, true));
  unexplored_nodes.push(make_shared<BPCNode>(*node, forbid_arcs_2, false));

  num_nodes_generated += 2;

  return true;
}

void BPCTree::print_result()
{
  instance->sylog.add_alg_result("stat", get_status_string());
  instance->sylog.add_alg_result("#node", BPCNode::num_explored);
  instance->sylog.add_alg_result("ub", ub);
  instance->sylog.add_alg_result("lb", lb);
  instance->sylog.add_alg_result("root_lb", root_lb);
  instance->sylog.add_alg_result("depth", max_depth);
  instance->sylog.add_alg_result("avg_time_sp", avg_time_on_sp);
  instance->sylog.add_alg_result("time", get_elapsed(tree_start_time));
  instance->sylog.add_alg_result("time_mp", total_time_on_mp);
  instance->sylog.add_alg_result("time_sp", total_time_on_sp);
  instance->sylog.add_alg_result("#spp", num_spp);
  instance->sylog.add_alg_result("time_sepa", total_time_on_sepa);
  instance->sylog.add_alg_result("#vc", num_vc_ctrs);
  instance->sylog.add_alg_result("#ac", num_ac_ctrs);
  instance->sylog.add_alg_result("#tc", num_tc_ctrs);
  instance->sylog.add_alg_result("#wc", num_wc_ctrs);
  instance->sylog.add_alg_result("#rc", num_rc_ctrs);
  instance->sylog.add_alg_result("#oc", num_oc_ctrs);
  instance->sylog.add_alg_result("#cc", num_cc_ctrs);
  instance->sylog.add_alg_result("#1c", num_1c_ctrs);
  instance->sylog.add_alg_result("#2c", num_2c_ctrs);
  instance->sylog.add_alg_result("#3c", num_3c_ctrs);
  instance->sylog.add_alg_result("#4c", num_4c_ctrs);
  instance->sylog.add_alg_result("time_heur", total_time_on_heur);
  instance->sylog.add_alg_result("time_br", total_time_on_branching);

  instance->sylog.print_result();
}

string BPCTree::get_status_string()
{
  switch (status)
  {
  case Status::ONGOING:
    return "ongoing";
    break;
  case Status::TIME_LIMIT_REACHED:
    return "time-limit-reached";
    break;
  case Status::OPTIMAL_FOUNDED:
    return "optimal";
    break;
  case Status::INFEASIBLE:
    return "infeasible";
    break;
  case Status::LP_OPTIMAL:
    return "lp_optimal";
    break;
  default:
    return "";
    break;
  }
}

void BPCTree::update_stats(shared_ptr<BPCNode> node)
{
  bool is_root = node->id == 1;

  // Update upper bound
  if (ub > node->ub)
  {
    ub = node->ub;
    node_attaining_ub = node;
    bound_updated = true;
  }

  // If is_root, then update lower bound and root gap
  if (is_root)
  {
    lb = node->lb;
    root_lb = lb;
    if (ub == DOUBLE_NULL)
    {
      gap = 1;
    }
    else
    {
      gap = (ub - lb) / ub;
    }
    gap_at_root = gap;
    bound_updated = false;
  }

  // Update gap if needed
  if (bound_updated)
  {
    if (ub == DOUBLE_NULL)
    {
      gap = 1;
      gap_at_root = 1;
    }
    else
    {
      gap = (ub - lb) / ub;
      gap_at_root = (ub - root_lb) / ub;
    }
    bound_updated = false;
  }

  // Update max_depth
  if (max_depth <= node->depth)
  {
    max_depth = node->depth;
  }

  // Update time records
  if (is_root)
  {
    avg_time_on_sp = node->avg_time_spent_on_sp;
    avg_time_on_erasing = node->time_on_erasing_arcs;
  }
  else
  {
    avg_time_on_sp = (total_time_on_sp + node->total_time_spent_on_sub) / ((total_time_on_sp / avg_time_on_sp) + node->all_times_spent_on_sub.size());

    avg_time_on_erasing = avg_time_on_erasing * (node->num_explored - 1) / (node->num_explored) + node->time_on_erasing_arcs / (node->num_explored);
  }
  elapsed_time = get_elapsed(tree_start_time);
  total_time_on_mp += node->total_time_spent_on_mas;
  total_time_on_sp += node->total_time_spent_on_sub;
  total_time_on_sepa += node->total_time_spent_on_sepa;
  total_time_on_heur += node->total_time_on_heur;

  num_vc_ctrs += node->num_vc_ctrs;
  num_ac_ctrs += node->num_ac_ctrs;
  num_tc_ctrs += node->num_tc_ctrs;
  num_wc_ctrs += node->num_wc_ctrs;
  num_rc_ctrs += node->num_rc_ctrs;
  num_oc_ctrs += node->num_oc_ctrs;
  num_cc_ctrs += node->num_cc_ctrs;
  num_1c_ctrs += node->num_1c_ctrs;
  num_2c_ctrs += node->num_2c_ctrs;
  num_3c_ctrs += node->num_3c_ctrs;
  num_4c_ctrs += node->num_4c_ctrs;

  num_spp += node->all_times_spent_on_sub.size();
}