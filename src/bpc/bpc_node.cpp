#include "bpc_node.h"
#include "sy_log.h"
#include "spp_teg.h"
#include "clique.h"

#include <chrono>
#include <numeric>
#include <memory>

int BPCNode::num_explored = 0;

/* Constructor that will be called to generate the root node */
BPCNode::BPCNode(
    shared_ptr<Instance> instance,
    shared_ptr<ColumnPool> pool,
    shared_ptr<MasterSolver> mas_solver,
    shared_ptr<SubSolver> sub_solver) : instance(instance), pool(pool), mas_solver(mas_solver), sub_solver(sub_solver)
{
  lb = DOUBLE_NULL;
  ub = DOUBLE_NULL;
  parent_lb = DOUBLE_NULL;
  depth = 0;
  id = 1;
  node_num_new_cols = 0;
  num_iter = 0;
  all_times_spent_on_sub = vector<double>();
  avg_time_spent_on_sp = 0.0;
  total_time_spent_on_sub = 0.0;
  total_time_spent_on_sepa = 0.0;
  total_time_spent_on_mas = 0.0;
  total_time_on_heur = 0.0;
  total_time_spent = 0.0;
  num_vc_ctrs = 0;
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

  // setup initial columns
  // one dummy
  auto dummy = make_shared<Column>(instance);
  dummy->make_dummy();
  pool->push_back(dummy);
  local_pool = *pool;
  cliques = make_shared<vector<Clique>>();

  // setup agent_order for primal heuristic
  agent_order = vector<tuple<int, int>>(); // id, |o - d|
  for (int agent_id = 0; agent_id < instance->param.num_agents; agent_id++)
  {

    agent_order.push_back(make_tuple(agent_id, instance->graph.get_vertex(Vertex::Type::DUMMY_ORIGIN, agent_id)->coord.dist_from(instance->graph.get_vertex(Vertex::Type::DUMMY_DESTINATION, agent_id)->coord)));
  }
  std::sort(agent_order.begin(), agent_order.end(), [](const tuple<int, int> &first, const tuple<int, int> &second)
            { return std::get<1>(first) < std::get<1>(second); }); 
}

/* Coonstructor that will be called from BPCTree::branch to generate child nodes */
BPCNode::BPCNode(const BPCNode &parent,
                 shared_ptr<BranchingRule> branching_rule,
                 bool is_left_child) : instance{parent.instance},
                                       pool{parent.pool},
                                       // local_pool{*pool},
                                       local_pool{parent.local_pool},
                                       parent_lb{parent.lb},
                                       erased_arcs{parent.erased_arcs},
                                       branching_rule{branching_rule},
                                       depth{parent.depth + 1},
                                       id{parent.id * 2 + (is_left_child ? 0 : 1)},
                                       mas_solver{parent.mas_solver},
                                       sub_solver{parent.sub_solver},
                                       cliques{parent.cliques},
                                       agent_order{parent.agent_order},
                                       ub{parent.ub}
{
  // initialize information using parent's information
  lb = DOUBLE_NULL;
  node_num_new_cols = 0;
  num_iter = 0;
  all_times_spent_on_sub = vector<double>();
  avg_time_spent_on_sp = 0.0;
  total_time_spent_on_sub = 0.0;
  total_time_spent_on_sepa = 0.0;
  total_time_spent_on_mas = 0.0;
  total_time_spent = 0.0;
  num_vc_ctrs = 0;
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
  update_erased_arcs();
  remove_incompatible_columns();
}

void BPCNode::update_erased_arcs()
{
  if (branching_rule == nullptr)
    return;
  branching_rule->add_erased_arcs(erased_arcs);
}

void BPCNode::remove_incompatible_columns()
{
  if (branching_rule == nullptr)
    return;

  ColumnPool new_local_pool;
  for (const auto &col : local_pool)
  {
    if (branching_rule->is_column_compatible(col))
    {
      new_local_pool.push_back(col);
    }
  }

  local_pool = new_local_pool;
}

/* Solve LP relaxation via column generation */
void BPCNode::solve()
{
  /* Initial setup */
  std::string ccg_status;
  auto node_start = get_now();
  mas_solver->reset_model();

  // remove trash value;
  base_columns = vector<tuple<shared_ptr<Column>, double>>();

  /* Solve master with dummy variables */
  auto mas_start = get_now();
  mas_result = mas_solver->solve_lp(local_pool);
  total_time_spent_on_mas += get_elapsed(mas_start);
#ifdef SHOW_MASSUB_RESULT
  println("            MP: obj {:.3f}  time {:.3f}", mas_result->obj_val, get_elapsed(node_start)); // 12 initial spaces
#endif

#ifdef COLGEN_DUAL_DUPLICATION_CHECK
  vector<MasterResult> mas_results;
  mas_results.push_back(mas_result);
#endif
  // iteration counter
  int num_cut_iter = 0;
  int num_col_iter = 0;

  /* Cut generation loop */
  int num_new_cut = 0;
  bool cut_added = true;
  while (cut_added)
  {
    if (get_elapsed(node_start) > instance->param.time_limit)
    {
      break;
    }
    num_col_iter = 0;
    /* Column generation loop */
    bool profitable_column_found = true;
    while (profitable_column_found)
    {
      if (get_elapsed(node_start) > instance->param.time_limit)
      {
        break;
      }

      /* Solve pricing subproblems */
      int num_new_cols = 0;
      if (!(num_cut_iter > 0 && num_col_iter == 0))
      {
        for (int agent_id = 0; agent_id < instance->param.num_agents; agent_id++)
        {
          auto sub_start = get_now();
          auto sub_result = sub_solver->solve(agent_id,
                                              local_pool,
                                              pool,
                                              mas_result,
                                              erased_arcs);
          auto sub_time = get_elapsed(sub_start);
          total_time_spent_on_sub += sub_time;
          all_times_spent_on_sub.push_back(sub_time);
          num_new_cols += sub_result.num_new_cols;
          node_num_new_cols += sub_result.num_new_cols;
        }
#ifdef SHOW_COLGEN_RESULT
        println("        Col-{}-{}: obj {:.3f}  |cols| {}  |new_cols| {}  time {:.3f}", num_cut_iter, num_col_iter, mas_result->obj_val, local_pool.size(), num_new_cols, get_elapsed(node_start)); // 8 initial spaces
#endif
      }

      if (num_new_cols > 0 || num_new_cut > 0) // some col or ctr added ==> re-solve master
      {
        num_iter += 1;
        num_col_iter += 1;
        /* Re-solve the Master with added columns */
        mas_start = get_now();
        mas_result = mas_solver->solve_lp(local_pool);
        total_time_spent_on_mas += get_elapsed(mas_start);

        num_new_cut = 0;

#ifdef COLGEN_DUAL_DUPLICATION_CHECK 
        for (auto result : mas_results)
        {
          if (std::equal(mas_result.cardi_duals.begin(), mas_result.cardi_duals.end(), result.cardi_duals.begin()))
          {
            println("same cardinality duals");
          }
        }
        mas_results.push_back(mas_result);
#endif
#ifdef SHOW_MASSUB_RESULT
        println("            MP: obj {:.3f}  time {:.3f}", mas_result->obj_val, get_elapsed(node_start));
        ; // 12 initial spaces
#endif
      }
      else // (num_new_cols <= 0) i.e. No profitable column exists
      {
        profitable_column_found = false;

        lb = mas_result->obj_val;
        base_columns = mas_result->solution;

        avg_time_spent_on_sp = std::accumulate(
                                   all_times_spent_on_sub.begin(),
                                   all_times_spent_on_sub.end(),
                                   0.0) /
                               all_times_spent_on_sub.size();
      }
    }

    if (instance->param.cut_depth >= 0 && depth > instance->param.cut_depth)
      break;

    if (instance->param.cut_rounds >= 0 && num_cut_iter > instance->param.cut_rounds)
      break;

    /* Check violated conflict constraints and add if exist */
    cut_added = false;
    if (get_elapsed(node_start) > instance->param.time_limit)
      continue;
    auto sepa_start = get_now();
    // Vertex conflict
    if (instance->param.ctrs_string.find('V') != std::string::npos)
    {
      // check violation by enumeration
      vt_values vertex_occu;
      sy_map<tuple<int, int>, vector<int>> vertex_occu_cols;
      int col_idx = 0;
      for (auto [col, val] : mas_result->solution)
      {
        int v_idx = 0;
        int time = 0;
        for (auto v : col->path)
        {
          if (v->type == Vertex::Type::NORMAL && v_idx <= col->path.size() - 2)
          {
            vertex_occu[make_tuple(v->id, time)] += val;
            vertex_occu_cols[make_tuple(v->id, time)].push_back(col_idx);
            time += 1;
          }
          v_idx += 1;
        }
        col_idx += 1;
      }
      // Add cuts if violated
      for (auto &[key, val] : vertex_occu)
      {
        if (val > 1 + CUT_VIOLATION_EPS)
        {
          cut_added = true;
          num_new_cut += 1;
          num_vc_ctrs += 1;
          for (auto col_idx : vertex_occu_cols[key])
          {
            auto [col, _] = mas_result->solution[col_idx];
            col->vertex_conflict_coeff[key] = 1;
            col->is_altered = true;
          }
        }
      }
    }
    // Arc conflict
    if (instance->param.ctrs_string.find('A') != std::string::npos)
    {
      // check violation by enumeration
      at_values arc_occu;
      sy_map<tuple<int, int>, vector<int>> arc_occu_cols;
      int col_idx = 0;
      for (auto [col, val] : mas_result->solution)
      {
        int v_idx = 0;
        int time = 0;
        for (auto v : col->path)
        {
          if (v->type == Vertex::Type::NORMAL && v_idx <= col->path.size() - 3)
          {
            auto w = col->path[v_idx + 1];
            if (v->id < w->id)
            {
              arc_occu[make_tuple(instance->graph.get_arc(v, w)->id, time)] += val;
              arc_occu_cols[make_tuple(instance->graph.get_arc(v, w)->id, time)].push_back(col_idx);
            }
            else if (v->id > w->id)
            {
              arc_occu[make_tuple(instance->graph.get_arc(w, v)->id, time)] += val;
              arc_occu_cols[make_tuple(instance->graph.get_arc(w, v)->id, time)].push_back(col_idx);
            }
            time += 1;
          }
          v_idx += 1;
        }
        col_idx += 1;
      }
      for (auto &[key, val] : arc_occu)
      {
        if (val > 1 + CUT_VIOLATION_EPS)
        {
          cut_added = true;
          num_new_cut += 1;
          num_ac_ctrs += 1;
          for (auto col_idx : arc_occu_cols[key])
          {
            auto [col, _] = mas_result->solution[col_idx];
            col->arc_conflict_coeff[key] = 1;
            col->is_altered = true;
          }
        }
      }
    }
    // Time conflict
    if (instance->param.ctrs_string.find('T') != std::string::npos)
    {
      // Check violation by enumeration
      ivtt_values time_occu;
      int col_idx = 0;
      for (auto [col, val] : mas_result->solution)
      {
        int time = 0;
        for (auto v : col->path)
        {
          if (v->type == Vertex::Type::NORMAL)
          {
            for (int delt = 0; delt <= instance->param.time_spacing; delt++)
            {
              time_occu[make_tuple(col->agent_id, v->id, time, delt)] += val;
            }
            for (int delt = 0; time - delt >= 0 && delt <= instance->param.time_spacing; delt++)
            {
              for (int j = 0; j < instance->param.num_agents; j++)
              {
                if (j != col->agent_id)
                {
                  time_occu[make_tuple(j, v->id, time - delt, delt)] += val;
                }
              }
            }
            time += 1;
          }
        }
        col_idx += 1;
      }
      // add cuts if violated
      for (auto &[key, val] : time_occu)
      {
        auto &[agent_id, vertex_id, time, delt] = key;
        if (val > 1 + CUT_VIOLATION_EPS)
        {
          mas_solver->time_conflict_is_needed[key] = true;
          cut_added = true;
          num_new_cut += 1;
          num_tc_ctrs += 1;
          for (auto col : *pool)
          {
            if (col->is_dummy)
              continue;
            bool is_in_cut = false;
            if (col->agent_id == agent_id && time + 1 <= col->path.size() - 1 && col->path[time + 1]->id == vertex_id)
            {
              is_in_cut = true;
            }
            else if (col->agent_id != agent_id)
            {
              if (time + delt + 1 <= col->path.size() - 1)
              {
                if (col->path[time + delt + 1]->id == vertex_id)
                  is_in_cut = true;
              }
            }
            if (is_in_cut)
            {
              col->time_conflict_coeff[key] = 1;
              col->is_altered = true;
            }
          }
        }
      }
    }
    // Wait conflict
    if (instance->param.ctrs_string.find('W') != std::string::npos)
    {
      // check violation by enumeration
      vt_values wait_occu;
      sy_map<tuple<int, int>, vector<int>> wait_occu_cols;
      int col_idx = 0;
      for (auto [col, val] : mas_result->solution)
      {
        int v_idx = 0;
        int time = 0;
        for (auto v : col->path)
        {
          if (v->type == Vertex::Type::NORMAL)
          {
            wait_occu[make_tuple(v->id, time)] += val;
            wait_occu_cols[make_tuple(v->id, time)].push_back(col_idx);
            if (v->id != col->path[v_idx - 1]->id)
            {
              wait_occu[make_tuple(v->id, time - 1)] += val;
              wait_occu_cols[make_tuple(v->id, time - 1)].push_back(col_idx);
            }
            time += 1;
          }
          v_idx += 1;
        }
        col_idx += 1;
      }
      // Add cuts if violated
      for (auto &[key, val] : wait_occu)
      {
        if (val > 1 + CUT_VIOLATION_EPS)
        {
          cut_added = true;
          num_new_cut += 1;
          num_wc_ctrs += 1;
          for (auto col_idx : wait_occu_cols[key])
          {
            auto [col, _] = mas_result->solution[col_idx];
            col->wait_conflict_coeff[key] = 1;
            col->is_altered = true;
          }
        }
      }
    }
    // Range conflict
    if (instance->param.ctrs_string.find('R') != std::string::npos)
    {
      // check violation by enumeration
      vt_values range_occu;
      sy_map<tuple<int, int>, vector<int>> range_occu_cols;
      int col_idx = 0;
      for (auto [col, val] : mas_result->solution)
      {
        int v_idx = 0;
        int time = 0;
        for (auto v : col->path)
        {
          if (v->type == Vertex::Type::NORMAL)
          {
            auto w = col->path[v_idx + 1];
            if (v->id != w->id)
            {
              for (int tau = 0; tau <= instance->param.time_spacing; tau++)
              {
                range_occu[make_tuple(v->id, time - tau)] += val;
                range_occu_cols[make_tuple(v->id, time - tau)].push_back(col_idx);
              }
            }
            else
            {
              range_occu[make_tuple(v->id, time - instance->param.time_spacing)] += val;
              range_occu_cols[make_tuple(v->id, time - instance->param.time_spacing)].push_back(col_idx);
            }
            time += 1;
          }
          v_idx += 1;
        }
        col_idx += 1;
      }
      // Add cuts if violated
      for (auto &[key, val] : range_occu)
      {
        if (val > 1 + CUT_VIOLATION_EPS)
        {
          cut_added = true;
          num_new_cut += 1;
          num_rc_ctrs += 1;
          sy_map<int, int> col_cnt;
          for (auto col_idx : range_occu_cols[key])
          {
            auto [col, _] = mas_result->solution[col_idx];
            col_cnt[col_idx] += 1;
            col->is_altered = true;
          }
          for (auto [col_idx, cnt] : col_cnt)
          {
            auto [col, _] = mas_result->solution[col_idx];
            col->range_conflict_coeff[key] = cnt;
          }
        }
      }
    }
    // Old time conflict
    if (instance->param.ctrs_string.find('O') != std::string::npos)
    {
      // Check violation by enumeration
      ivt_values old_occu;
      int col_idx = 0;
      for (auto [col, val] : mas_result->solution)
      {
        int time = 0;
        for (auto v : col->path)
        {
          if (v->type == Vertex::Type::NORMAL)
          {
            old_occu[make_tuple(col->agent_id, v->id, time)] += val;
            for (int delt = 0; time - delt >= 0 && delt <= instance->param.time_spacing; delt++)
            {
              for (int j = 0; j < instance->param.num_agents; j++)
              {
                if (j != col->agent_id)
                {
                  old_occu[make_tuple(j, v->id, time - delt)] += val / ((double)instance->param.time_spacing + 1);
                }
              }
            }
            time += 1;
          }
        }
        col_idx += 1;
      }
      // add cuts if violated
      for (auto &[key, val] : old_occu)
      {
        auto &[agent_id, vertex_id, time] = key;
        if (val > 1 + CUT_VIOLATION_EPS)
        {
          mas_solver->old_conflict_is_needed[key] = true;
          cut_added = true;
          num_new_cut += 1;
          num_oc_ctrs += 1;
          
          for (auto col : *pool)
          {
            if (col->is_dummy)
              continue;
            bool is_in_cut = false;
            double coeff = 0;
            if (col->agent_id == agent_id && time + 1 <= col->path.size() - 1 && col->path[time + 1]->id == vertex_id)
            {
              is_in_cut = true;
              coeff += 1;
            }
            else if (col->agent_id != agent_id)
            {
              for (int delt = 0; delt <= instance->param.time_spacing; delt++)
              {
                if (time + delt + 1 <= col->path.size() - 1)
                {
                  if (col->path[time + delt + 1]->id == vertex_id)
                  {
                    is_in_cut = true;
                    coeff += 1 / ((double)instance->param.time_spacing + 1);
                  }
                }
              }
            }
            if (is_in_cut)
            {
              col->old_conflict_coeff[key] = coeff;
              col->is_altered = true;
            }
          }
        }
      }
    }
    // All TVC conflicts
    if (instance->param.ctrs_string.find("1234") != std::string::npos)
    {
      // Check violation by enumeration
      ivwt_values tvc1_occu, tvc2_occu, tvc3_occu, tvc4_occu;

      // option 1: check violated i,v,w,t for candidate i,v,w
      //// inspect i,v,w
      sy_map<tuple<int, int, int>, double> tvc_apx_occu;
      for (auto [col, val] : mas_result->solution)
      {
        int v_idx = -1;
        for (auto v : col->path)
        {
          v_idx += 1;
          if (v->type == Vertex::Type::NORMAL)
          {
            auto w = col->path[v_idx + 1];
            if (v == w)
            {
              for (auto a : instance->graph.get_outout_arcs(v))
              {
                auto u = a->to;
                if (u->type == Vertex::Type::NORMAL)
                {
                  tvc_apx_occu[make_tuple(col->agent_id, v->id, u->id)] += val;
                  tvc_apx_occu[make_tuple(col->agent_id, u->id, v->id)] += val;
                }
              }
              continue;
            }
            else if (w->type == Vertex::Type::NORMAL) // v != w
            {
              // i == col->agent_id
              tvc_apx_occu[make_tuple(col->agent_id, v->id, w->id)] += val;
              tvc_apx_occu[make_tuple(col->agent_id, w->id, v->id)] += val;
              for (auto a : instance->graph.get_outout_arcs(v))
              {
                auto u = a->to;
                if (u == w)
                  continue;
                tvc_apx_occu[make_tuple(col->agent_id, v->id, u->id)] += val;
                tvc_apx_occu[make_tuple(col->agent_id, u->id, v->id)] += val;
              }
              for (auto a : instance->graph.get_outout_arcs(w))
              {
                auto u = a->to;
                if (u == v)
                  continue;
                tvc_apx_occu[make_tuple(col->agent_id, w->id, u->id)] += val;
                tvc_apx_occu[make_tuple(col->agent_id, u->id, w->id)] += val;
              }
              // i != col->agent_id
              for (int j = 0; j < instance->param.num_agents; j++)
              {
                if (j != col->agent_id)
                {
                  tvc_apx_occu[make_tuple(j, v->id, w->id)] += val;
                  tvc_apx_occu[make_tuple(j, w->id, v->id)] += val;
                }
              }
            }
          }
        }
      }
      //// exact evaluation for candidate i,v,w for each t
      for (auto [col, val] : mas_result->solution)
      {
        int time = 0;
        int v_idx = -1;
        for (auto v : col->path)
        {
          v_idx += 1;
          if (v->type == Vertex::Type::NORMAL)
          {
            auto w = col->path[v_idx + 1];
            if (v == w)
            {
              for (auto a : instance->graph.get_outout_arcs(v))
              {
                auto u = a->to;
                if (u->type == Vertex::Type::NORMAL && tvc_apx_occu[make_tuple(col->agent_id, v->id, u->id)] > 1 + CUT_VIOLATION_EPS)
                {
                  tvc1_occu[make_tuple(col->agent_id, std::max(v->id, u->id), std::min(v->id, u->id), time)] += val;
                  tvc2_occu[make_tuple(col->agent_id, std::max(v->id, u->id), std::min(v->id, u->id), time)] += val;
                  tvc3_occu[make_tuple(col->agent_id, v->id, u->id, time)] += val;
                  tvc3_occu[make_tuple(col->agent_id, u->id, v->id, time)] += val;
                  tvc4_occu[make_tuple(col->agent_id, v->id, u->id, time)] += val;
                  tvc4_occu[make_tuple(col->agent_id, u->id, v->id, time)] += val;
                }
              }
              time += 1;
              continue;
            }
            else if (w->type == Vertex::Type::NORMAL)
            {
              // i == col->agent_id
              if (tvc_apx_occu[make_tuple(col->agent_id, v->id, w->id)] > 1 + CUT_VIOLATION_EPS)
              {
                tvc1_occu[make_tuple(col->agent_id, std::max(v->id, w->id), std::min(v->id, w->id), time)] += val;
                tvc2_occu[make_tuple(col->agent_id, std::max(v->id, w->id), std::min(v->id, w->id), time)] += val;
                tvc3_occu[make_tuple(col->agent_id, v->id, w->id, time)] += val;
                tvc3_occu[make_tuple(col->agent_id, w->id, v->id, time)] += val;
                tvc4_occu[make_tuple(col->agent_id, v->id, w->id, time)] += val;
                tvc4_occu[make_tuple(col->agent_id, w->id, v->id, time)] += val;
              }
              for (auto a : instance->graph.get_outout_arcs(v))
              {
                auto u = a->to;
                if (u != w && tvc_apx_occu[make_tuple(col->agent_id, v->id, u->id)] > 1 + CUT_VIOLATION_EPS)
                {
                  tvc1_occu[make_tuple(col->agent_id, std::max(v->id, u->id), std::min(v->id, u->id), time)] += val;
                  tvc2_occu[make_tuple(col->agent_id, std::max(v->id, u->id), std::min(v->id, u->id), time)] += val;
                  tvc3_occu[make_tuple(col->agent_id, v->id, u->id, time)] += val;
                  tvc4_occu[make_tuple(col->agent_id, v->id, u->id, time)] += val;
                }
              }
              for (auto a : instance->graph.get_outout_arcs(w))
              {
                auto u = a->to;
                if (u != v && tvc_apx_occu[make_tuple(col->agent_id, w->id, u->id)] > 1 + CUT_VIOLATION_EPS)
                {
                  tvc1_occu[make_tuple(col->agent_id, std::max(w->id, u->id), std::min(w->id, u->id), time)] += val;
                  tvc2_occu[make_tuple(col->agent_id, std::max(w->id, u->id), std::min(w->id, u->id), time)] += val;
                  tvc3_occu[make_tuple(col->agent_id, w->id, u->id, time)] += val;
                }
              }
              for (auto a : instance->graph.get_inin_arcs(v))
              {
                auto u = a->from;
                if (u != w && tvc_apx_occu[make_tuple(col->agent_id, v->id, u->id)] > 1 + CUT_VIOLATION_EPS)
                  tvc4_occu[make_tuple(col->agent_id, u->id, v->id, time)] += val;
              }
              for (auto a : instance->graph.get_inin_arcs(w))
              {
                auto u = a->from;
                if (u != v && tvc_apx_occu[make_tuple(col->agent_id, w->id, u->id)] > 1 + CUT_VIOLATION_EPS)
                {
                  tvc3_occu[make_tuple(col->agent_id, u->id, w->id, time)] += val;
                  tvc4_occu[make_tuple(col->agent_id, u->id, w->id, time)] += val;
                }
              }
              // i != col->agent_id
              for (int j = 0; j < instance->param.num_agents; j++)
              {
                if (j != col->agent_id && tvc_apx_occu[make_tuple(j, v->id, w->id)] > 1 + CUT_VIOLATION_EPS)
                {
                  tvc1_occu[make_tuple(j, std::max(v->id, w->id), std::min(v->id, w->id), time)] += val;
                  tvc2_occu[make_tuple(j, std::max(v->id, w->id), std::min(v->id, w->id), time)] += val;
                  tvc3_occu[make_tuple(j, v->id, w->id, time)] += val;
                  tvc4_occu[make_tuple(j, v->id, w->id, time)] += val;
                  for (int delt = 1; delt <= instance->param.time_spacing - 1; delt++)
                  {
                    if (time - delt >= 0)
                    {
                      tvc1_occu[make_tuple(j, std::max(v->id, w->id), std::min(v->id, w->id), time - delt)] += val;
                      tvc3_occu[make_tuple(j, v->id, w->id, time - delt)] += val;
                      tvc3_occu[make_tuple(j, w->id, v->id, time - delt)] += val;
                    }
                    tvc2_occu[make_tuple(j, std::max(v->id, w->id), std::min(v->id, w->id), time + delt)] += val;
                    tvc4_occu[make_tuple(j, v->id, w->id, time + delt)] += val;
                    tvc4_occu[make_tuple(j, w->id, v->id, time + delt)] += val;
                  }
                  if (time - instance->param.time_spacing >= 0)
                    tvc3_occu[make_tuple(j, v->id, w->id, time - instance->param.time_spacing)] += val;
                  tvc4_occu[make_tuple(j, v->id, w->id, time + instance->param.time_spacing)] += val;
                }
              }
            }
            time += 1;
          }
        }
      }

      // add cuts if violated
      for (auto &[key, val] : tvc1_occu)
      {
        auto &[agent_id, v_id, w_id, time] = key;
        if (v_id <= w_id)
          continue;
        if (val > 1 + CUT_VIOLATION_EPS)
        {
          mas_solver->tvc1_conflict_is_needed[key] = true;
          cut_added = true;
          num_new_cut += 1;
          num_1c_ctrs += 1;

          for (auto col : *pool)
          {
            if (col->is_dummy)
              continue;
            double coeff = 0;
            if (col->agent_id == agent_id && time + 1 + 1 <= col->path.size() - 1)
            {
              if (col->path[time + 1]->id == v_id || col->path[time + 1]->id == w_id || col->path[time + 1 + 1]->id == v_id || col->path[time + 1 + 1]->id == w_id)
              {
                coeff += 1;
              }
            }
            else if (col->agent_id != agent_id)
            {
              for (int delt = 0; delt <= instance->param.time_spacing - 1; delt++)
              {
                if (time + 1 + delt + 1 <= col->path.size() - 1)
                {
                  if ((col->path[time + 1 + delt]->id == v_id && col->path[time + 1 + delt + 1]->id == w_id) || (col->path[time + 1 + delt]->id == w_id && col->path[time + 1 + delt + 1]->id == v_id))
                  {
                    coeff += 1;
                  }
                }
              }
            }
            if (coeff > 0)
            {
              col->tvc1_conflict_coeff[key] = coeff;
              col->is_altered = true;
            }
          }
        }
      }
      for (auto &[key, val] : tvc2_occu)
      {
        auto &[agent_id, v_id, w_id, time] = key;
        if (v_id <= w_id)
          continue;
        if (val > 1 + CUT_VIOLATION_EPS)
        {
          mas_solver->tvc2_conflict_is_needed[key] = true;
          cut_added = true;
          num_new_cut += 1;
          num_2c_ctrs += 1;

          for (auto col : *pool)
          {
            if (col->is_dummy)
              continue;
            double coeff = 0;
            if (col->agent_id == agent_id && time + 1 + 1 <= col->path.size() - 1)
            {
              if (col->path[time + 1]->id == v_id || col->path[time + 1]->id == w_id || col->path[time + 1 + 1]->id == v_id || col->path[time + 1 + 1]->id == w_id)
              {
                coeff += 1;
              }
            }
            else if (col->agent_id != agent_id)
            {
              for (int delt = 0; delt <= instance->param.time_spacing - 1; delt++)
              {
                if (time + 1 - delt >= 0 && time + 1 - delt + 1 <= col->path.size() - 1)
                {
                  if ((col->path[time + 1 - delt]->id == v_id && col->path[time + 1 - delt + 1]->id == w_id) || (col->path[time + 1 - delt]->id == w_id && col->path[time + 1 - delt + 1]->id == v_id))
                  {
                    coeff += 1;
                  }
                }
              }
            }
            if (coeff > 0)
            {
              col->tvc2_conflict_coeff[key] = coeff;
              col->is_altered = true;
            }
          }
        }
      }
      for (auto &[key, val] : tvc3_occu)
      {
        auto &[agent_id, v_id, w_id, time] = key;

        if (val > 1 + CUT_VIOLATION_EPS)
        {
          mas_solver->tvc3_conflict_is_needed[key] = true;
          cut_added = true;
          num_new_cut += 1;
          num_3c_ctrs += 1;
          for (auto col : *pool)
          {
            if (col->is_dummy)
              continue;
            double coeff = 0;
            if (col->agent_id == agent_id && time + 1 + 1 <= col->path.size() - 1)
            {
              if (col->path[time + 1]->id == v_id || col->path[time + 2]->id == v_id || col->path[time + 2]->id == w_id)
              {
                coeff += 1;
              }
            }
            else if (col->agent_id != agent_id)
            {
              if (time + 1 + 1 <= col->path.size() - 1)
                if ((col->path[time + 1]->id == v_id && col->path[time + 2]->id == w_id))
                {
                  coeff += 1;
                }
              for (int delt = 1; delt <= instance->param.time_spacing - 1; delt++)
              {
                if (time + 1 + delt + 1 <= col->path.size() - 1)
                {
                  if ((col->path[time + 1 + delt]->id == v_id && col->path[time + 1 + delt + 1]->id == w_id) || (col->path[time + 1 + delt]->id == w_id && col->path[time + 1 + delt + 1]->id == v_id))
                  {
                    coeff += 1;
                  }
                }
              }
              if (time + 2 + instance->param.time_spacing <= col->path.size() - 1)
              {
                if ((col->path[time + 1 + instance->param.time_spacing]->id == v_id && col->path[time + 1 + instance->param.time_spacing + 1]->id == w_id))
                {
                  coeff += 1;
                }
              }
            }
            if (coeff > 0)
            {
              col->tvc3_conflict_coeff[key] = coeff;
              col->is_altered = true;
            }
          }
        }
      }
      for (auto &[key, val] : tvc4_occu)
      {
        auto &[agent_id, v_id, w_id, time] = key;

        if (val > 1 + CUT_VIOLATION_EPS)
        {
          mas_solver->tvc4_conflict_is_needed[key] = true;
          cut_added = true;
          num_new_cut += 1;
          num_4c_ctrs += 1;

          for (auto col : *pool)
          {
            if (col->is_dummy)
              continue;
            double coeff = 0;
            if (col->agent_id == agent_id && time + 1 + 1 <= col->path.size() - 1)
            {
              if (col->path[time + 1]->id == v_id || col->path[time + 1]->id == w_id || col->path[time + 2]->id == w_id)
              {
                coeff += 1;
              }
            }
            else if (col->agent_id != agent_id)
            {
              if (time + 1 + 1 <= col->path.size() - 1)
                if ((col->path[time + 1]->id == v_id && col->path[time + 2]->id == w_id))
                {
                  coeff += 1;
                }
              for (int delt = 1; delt <= instance->param.time_spacing - 1; delt++)
              {
                if (time + 1 - delt + 1 <= col->path.size() - 1 && time + 1 - delt >= 0)
                {
                  if ((col->path[time + 1 - delt]->id == v_id && col->path[time + 1 - delt + 1]->id == w_id) || (col->path[time + 1 - delt]->id == w_id && col->path[time + 1 - delt + 1]->id == v_id))
                  {
                    coeff += 1;
                  }
                }
              }
              if (time + 2 - instance->param.time_spacing <= col->path.size() - 1 && time + 1 - instance->param.time_spacing >= 0)
              {
                if ((col->path[time + 1 - instance->param.time_spacing]->id == v_id && col->path[time + 1 - instance->param.time_spacing + 1]->id == w_id))
                {
                  coeff += 1;
                }
              }
            }
            if (coeff > 0)
            {
              col->tvc4_conflict_coeff[key] = coeff;
              col->is_altered = true;
            }
          }
        }
      }
    }
    else
    {
      // 1c / Two-vertex-clique-1 conflict
      if (instance->param.ctrs_string.find('1') != std::string::npos)
      {
        // Check violation by enumeration
        ivwt_values tvc1_occu;
        int col_idx = 0;
        for (auto [col, val] : mas_result->solution)
        {
          int time = 0;
          int v_idx = -1;
          for (auto v : col->path)
          {
            v_idx += 1;
            if (v->type == Vertex::Type::NORMAL)
            {
              auto w = col->path[v_idx + 1];
              if (v == w)
              {
                for (auto a : instance->graph.get_outout_arcs(v))
                {
                  if (a->to->type == Vertex::Type::NORMAL)
                  {
                    auto u = a->to;
                    tvc1_occu[make_tuple(col->agent_id, v->id, u->id, time)] += val;
                    tvc1_occu[make_tuple(col->agent_id, u->id, v->id, time)] += val;
                  }
                }
                time += 1;
                continue;
              }
              else if (v != w)
              {
                for (auto a : instance->graph.get_outout_arcs(v))
                {
                  auto u = a->to;
                  tvc1_occu[make_tuple(col->agent_id, v->id, u->id, time)] += val;
                }
                for (auto a : instance->graph.get_outout_arcs(w))
                {
                  auto u = a->to;
                  tvc1_occu[make_tuple(col->agent_id, w->id, u->id, time)] += val;
                }
                for (auto a : instance->graph.get_inin_arcs(v))
                {
                  if (a->from == v || a->from == w)
                    continue;
                  auto u = a->from;
                  tvc1_occu[make_tuple(col->agent_id, u->id, v->id, time)] += val;
                }
                for (auto a : instance->graph.get_inin_arcs(w))
                {
                  if (a->from == v || a->from == w)
                    continue;
                  auto u = a->from;
                  tvc1_occu[make_tuple(col->agent_id, u->id, w->id, time)] += val;
                }
                for (int j = 0; j < instance->param.num_agents; j++)
                {
                  if (j != col->agent_id)
                  {
                    for (int delt = 0; time - delt >= 0 && delt <= instance->param.time_spacing - 1; delt++)
                    {
                      tvc1_occu[make_tuple(j, v->id, w->id, time - delt)] += val;
                      tvc1_occu[make_tuple(j, w->id, v->id, time - delt)] += val;
                    }
                  }
                }
              }
              time += 1;
            }
          }
          col_idx += 1;
        }
        // add cuts if violated
        for (auto &[key, val] : tvc1_occu)
        {
          auto &[agent_id, v_id, w_id, time] = key;
          if (v_id <= w_id)
            continue;
          if (val > 1 + CUT_VIOLATION_EPS)
          {
            mas_solver->tvc1_conflict_is_needed[key] = true;
            cut_added = true;
            num_new_cut += 1;
            num_1c_ctrs += 1;
            for (auto col : *pool)
            {
              if (col->is_dummy)
                continue;
              bool is_in_cut = false;
              double coeff = 0;
              if (col->agent_id == agent_id && time + 1 + 1 <= col->path.size() - 1)
              {
                if (col->path[time + 1]->id == v_id || col->path[time + 1]->id == w_id || col->path[time + 1 + 1]->id == v_id || col->path[time + 1 + 1]->id == w_id)
                {
                  is_in_cut = true;
                  coeff += 1;
                }
              }
              else if (col->agent_id != agent_id)
              {
                for (int delt = 0; delt <= instance->param.time_spacing - 1; delt++)
                {
                  if (time + 1 + delt + 1 <= col->path.size() - 1)
                  {
                    if ((col->path[time + 1 + delt]->id == v_id && col->path[time + 1 + delt + 1]->id == w_id) || (col->path[time + 1 + delt]->id == w_id && col->path[time + 1 + delt + 1]->id == v_id))
                    {
                      is_in_cut = true;
                      coeff += 1;
                    }
                  }
                }
              }
              if (is_in_cut)
              {
                col->tvc1_conflict_coeff[key] = coeff;
                col->is_altered = true;
              }
            }
          }
        }
      }
      // 2c / Two-vertex-clique-2 conflict
      if (instance->param.ctrs_string.find('2') != std::string::npos)
      {
        // Check violation by enumeration
        ivwt_values tvc2_occu;
        int col_idx = 0;
        for (auto [col, val] : mas_result->solution)
        {
          int time = 0;
          int v_idx = -1;
          for (auto v : col->path)
          {
            v_idx += 1;
            if (v->type == Vertex::Type::NORMAL)
            {
              auto w = col->path[v_idx + 1];
              if (v == w)
              {
                for (auto a : instance->graph.get_outout_arcs(v))
                {
                  if (a->to->type == Vertex::Type::NORMAL)
                  {
                    auto u = a->to;
                    tvc2_occu[make_tuple(col->agent_id, v->id, u->id, time)] += val;
                    tvc2_occu[make_tuple(col->agent_id, u->id, v->id, time)] += val;
                  }
                }
                time += 1;
                continue;
              }
              else if (v != w)
              {
                for (auto a : instance->graph.get_outout_arcs(v))
                {
                  auto u = a->to;
                  tvc2_occu[make_tuple(col->agent_id, v->id, u->id, time)] += val;
                }
                for (auto a : instance->graph.get_outout_arcs(w))
                {
                  auto u = a->to;
                  tvc2_occu[make_tuple(col->agent_id, w->id, u->id, time)] += val;
                }
                for (auto a : instance->graph.get_inin_arcs(v))
                {
                  if (a->from == v || a->from == w)
                    continue;
                  auto u = a->from;
                  tvc2_occu[make_tuple(col->agent_id, u->id, v->id, time)] += val;
                }
                for (auto a : instance->graph.get_inin_arcs(w))
                {
                  if (a->from == v || a->from == w)
                    continue;
                  auto u = a->from;
                  tvc2_occu[make_tuple(col->agent_id, u->id, w->id, time)] += val;
                }
                for (int j = 0; j < instance->param.num_agents; j++)
                {
                  if (j != col->agent_id)
                  {
                    for (int delt = 0; delt <= instance->param.time_spacing - 1; delt++)
                    {
                      tvc2_occu[make_tuple(j, v->id, w->id, time + delt)] += val;
                      tvc2_occu[make_tuple(j, w->id, v->id, time + delt)] += val;
                    }
                  }
                }
              }
              time += 1;
            }
          }
          col_idx += 1;
        }
        // add cuts if violated
        for (auto &[key, val] : tvc2_occu)
        {
          auto &[agent_id, v_id, w_id, time] = key;
          if (v_id <= w_id)
            continue;
          if (val > 1 + CUT_VIOLATION_EPS)
          {
            mas_solver->tvc2_conflict_is_needed[key] = true;
            cut_added = true;
            num_new_cut += 1;
            num_2c_ctrs += 1;

            for (auto col : *pool)
            {
              if (col->is_dummy)
                continue;
              bool is_in_cut = false;
              double coeff = 0;
              if (col->agent_id == agent_id && time + 1 + 1 <= col->path.size() - 1)
              {
                if (col->path[time + 1]->id == v_id || col->path[time + 1]->id == w_id || col->path[time + 1 + 1]->id == v_id || col->path[time + 1 + 1]->id == w_id)
                {
                  is_in_cut = true;
                  coeff += 1;
                }
              }
              else if (col->agent_id != agent_id)
              {
                for (int delt = 0; delt <= instance->param.time_spacing - 1; delt++)
                {
                  if (time + 1 - delt >= 0 && time + 1 - delt + 1 <= col->path.size() - 1)
                  {
                    if ((col->path[time + 1 - delt]->id == v_id && col->path[time + 1 - delt + 1]->id == w_id) || (col->path[time + 1 - delt]->id == w_id && col->path[time + 1 - delt + 1]->id == v_id))
                    {
                      is_in_cut = true;
                      coeff += 1;
                    }
                  }
                }
              }
              if (is_in_cut)
              {
                col->tvc2_conflict_coeff[key] = coeff;
                col->is_altered = true;
              }
            }
          }
        }
      }
      // 3c / Two-vertex-clique-3 conflict
      if (instance->param.ctrs_string.find('3') != std::string::npos)
      {
        // Check violation by enumeration
        ivwt_values tvc3_occu;
        int col_idx = 0;
        for (auto [col, val] : mas_result->solution)
        {
          int time = 0;
          int v_idx = -1;
          for (auto v : col->path)
          {
            v_idx += 1;
            if (v->type == Vertex::Type::NORMAL)
            {
              auto w = col->path[v_idx + 1];
              if (v == w)
              {
                for (auto a : instance->graph.get_outout_arcs(v))
                {
                  if (a->to->type == Vertex::Type::NORMAL)
                  {
                    auto u = a->to;
                    tvc3_occu[make_tuple(col->agent_id, v->id, u->id, time)] += val;
                    tvc3_occu[make_tuple(col->agent_id, u->id, v->id, time)] += val;
                  }
                }
                time += 1;
                continue;
              }
              else if (v != w)
              {
                for (auto a : instance->graph.get_outout_arcs(v))
                {
                  auto u = a->to;
                  tvc3_occu[make_tuple(col->agent_id, v->id, u->id, time)] += val;
                }
                for (auto a : instance->graph.get_outout_arcs(w))
                {
                  auto u = a->to;
                  tvc3_occu[make_tuple(col->agent_id, w->id, u->id, time)] += val;
                }
                for (auto a : instance->graph.get_inin_arcs(w))
                {
                  if (a->from == v || a->from == w)
                    continue;
                  auto u = a->from;
                  tvc3_occu[make_tuple(col->agent_id, u->id, w->id, time)] += val;
                }
                for (int j = 0; j < instance->param.num_agents; j++)
                {
                  if (j != col->agent_id)
                  {
                    tvc3_occu[make_tuple(j, v->id, w->id, time)] += val;
                    for (int delt = 1; time - delt >= 0 && delt <= instance->param.time_spacing - 1; delt++)
                    {
                      tvc3_occu[make_tuple(j, v->id, w->id, time - delt)] += val;
                      tvc3_occu[make_tuple(j, w->id, v->id, time - delt)] += val;
                    }
                    tvc3_occu[make_tuple(j, v->id, w->id, time - instance->param.time_spacing)] += val;
                  }
                }
              }
              time += 1;
            }
          }
          col_idx += 1;
        }
        // add cuts if violated
        for (auto &[key, val] : tvc3_occu)
        {
          auto &[agent_id, v_id, w_id, time] = key;

          if (val > 1 + CUT_VIOLATION_EPS)
          {
            mas_solver->tvc3_conflict_is_needed[key] = true;
            cut_added = true;
            num_new_cut += 1;
            num_3c_ctrs += 1;

            for (auto col : *pool)
            {
              if (col->is_dummy)
                continue;
              bool is_in_cut = false;
              double coeff = 0;
              if (col->agent_id == agent_id && time + 1 + 1 <= col->path.size() - 1)
              {
                if (col->path[time + 1]->id == v_id || col->path[time + 2]->id == v_id || col->path[time + 2]->id == w_id)
                {
                  is_in_cut = true;
                  coeff += 1;
                }
              }
              else if (col->agent_id != agent_id)
              {
                if (time + 1 + 1 <= col->path.size() - 1)
                  if ((col->path[time + 1]->id == v_id && col->path[time + 2]->id == w_id))
                  {
                    is_in_cut = true;
                    coeff += 1;
                  }
                for (int delt = 1; delt <= instance->param.time_spacing - 1; delt++)
                {
                  if (time + 1 + delt + 1 <= col->path.size() - 1)
                  {
                    if ((col->path[time + 1 + delt]->id == v_id && col->path[time + 1 + delt + 1]->id == w_id) || (col->path[time + 1 + delt]->id == w_id && col->path[time + 1 + delt + 1]->id == v_id))
                    {
                      is_in_cut = true;
                      coeff += 1;
                    }
                  }
                }
                if (time + 2 + instance->param.time_spacing <= col->path.size() - 1)
                {
                  if ((col->path[time + 1 + instance->param.time_spacing]->id == v_id && col->path[time + 1 + instance->param.time_spacing + 1]->id == w_id))
                  {
                    is_in_cut = true;
                    coeff += 1;
                  }
                }
              }
              if (is_in_cut)
              {
                col->tvc3_conflict_coeff[key] = coeff;
                col->is_altered = true;
              }
            }
          }
        }
      }
      // 4c / Two-vertex-clique-4 conflict
      if (instance->param.ctrs_string.find('4') != std::string::npos)
      {
        // Check violation by enumeration
        ivwt_values tvc4_occu;
        int col_idx = 0;
        for (auto [col, val] : mas_result->solution)
        {
          int time = 0;
          int v_idx = -1;
          for (auto v : col->path)
          {
            v_idx += 1;
            if (v->type == Vertex::Type::NORMAL)
            {
              auto w = col->path[v_idx + 1];
              if (v == w)
              {
                for (auto a : instance->graph.get_outout_arcs(v))
                {
                  if (a->to->type == Vertex::Type::NORMAL)
                  {
                    auto u = a->to;
                    tvc4_occu[make_tuple(col->agent_id, v->id, u->id, time)] += val;
                    tvc4_occu[make_tuple(col->agent_id, u->id, v->id, time)] += val;
                  }
                }
                time += 1;
                continue;
              }
              else if (v != w)
              {
                for (auto a : instance->graph.get_outout_arcs(v))
                {
                  auto u = a->to;
                  tvc4_occu[make_tuple(col->agent_id, v->id, u->id, time)] += val;
                }
                for (auto a : instance->graph.get_outout_arcs(w))
                {
                  auto u = a->to;
                  if (u == v)
                  {
                    tvc4_occu[make_tuple(col->agent_id, w->id, u->id, time)] += val;
                  }
                }
                for (auto a : instance->graph.get_inin_arcs(v))
                {
                  if (a->from == v || a->from == w)
                    continue;
                  auto u = a->from;
                  tvc4_occu[make_tuple(col->agent_id, u->id, v->id, time)] += val;
                }
                for (auto a : instance->graph.get_inin_arcs(w))
                {
                  if (a->from == v || a->from == w)
                    continue;
                  auto u = a->from;
                  tvc4_occu[make_tuple(col->agent_id, u->id, w->id, time)] += val;
                }
                for (int j = 0; j < instance->param.num_agents; j++)
                {
                  if (j != col->agent_id)
                  {
                    tvc4_occu[make_tuple(j, v->id, w->id, time)] += val;
                    for (int delt = 1; delt <= instance->param.time_spacing - 1; delt++)
                    {
                      tvc4_occu[make_tuple(j, v->id, w->id, time + delt)] += val;
                      tvc4_occu[make_tuple(j, w->id, v->id, time + delt)] += val;
                    }
                    tvc4_occu[make_tuple(j, v->id, w->id, time + instance->param.time_spacing)] += val;
                  }
                }
              }
              time += 1;
            }
          }
          col_idx += 1;
        }
        // add cuts if violated
        for (auto &[key, val] : tvc4_occu)
        {
          auto &[agent_id, v_id, w_id, time] = key;

          if (val > 1 + CUT_VIOLATION_EPS)
          {
            mas_solver->tvc4_conflict_is_needed[key] = true;
            cut_added = true;
            num_new_cut += 1;
            num_4c_ctrs += 1;
            
            for (auto col : *pool)
            {
              if (col->is_dummy)
                continue;
              bool is_in_cut = false;
              double coeff = 0;
              if (col->agent_id == agent_id && time + 1 + 1 <= col->path.size() - 1)
              {
                if (col->path[time + 1]->id == v_id || col->path[time + 1]->id == w_id || col->path[time + 2]->id == w_id)
                {
                  is_in_cut = true;
                  coeff += 1;
                }
              }
              else if (col->agent_id != agent_id)
              {
                if (time + 1 + 1 <= col->path.size() - 1)
                  if ((col->path[time + 1]->id == v_id && col->path[time + 2]->id == w_id))
                  {
                    is_in_cut = true;
                    coeff += 1;
                  }
                for (int delt = 1; delt <= instance->param.time_spacing - 1; delt++)
                {
                  if (time + 1 - delt + 1 <= col->path.size() - 1 && time + 1 - delt >= 0)
                  {
                    if ((col->path[time + 1 - delt]->id == v_id && col->path[time + 1 - delt + 1]->id == w_id) || (col->path[time + 1 - delt]->id == w_id && col->path[time + 1 - delt + 1]->id == v_id))
                    {
                      is_in_cut = true;
                      coeff += 1;
                    }
                  }
                }
                if (time + 2 - instance->param.time_spacing <= col->path.size() - 1 && time + 1 - instance->param.time_spacing >= 0)
                {
                  if ((col->path[time + 1 - instance->param.time_spacing]->id == v_id && col->path[time + 1 - instance->param.time_spacing + 1]->id == w_id))
                  {
                    is_in_cut = true;
                    coeff += 1;
                  }
                }
              }
              if (is_in_cut)
              {
                col->tvc4_conflict_coeff[key] = coeff;
                col->is_altered = true;
                // println("  c2 {} : {}", col_idx, col->name());
              }
            }
          }
        }
      }
    }
    // Clique conflict
    if (instance->param.ctrs_string.find('C') != std::string::npos && depth <= 0)
    {
      /* Initialize out clique_options */
      clique_options *opts = static_cast<clique_options *>(malloc(sizeof(clique_options)));
      // opts->time_function = NULL;
      opts->time_function = clique_print_time_with_limit;
      opts->reorder_function = NULL;
      opts->reorder_map = NULL;
      opts->output = stdout;
      opts->user_data = NULL;
      opts->user_function = NULL;
      opts->clique_list_length = 100000;
      opts->clique_list = static_cast<set_t *>(malloc(sizeof(set_t) * opts->clique_list_length));

      // 0. generate incompatibility graph
      int cl_threshold = 1000;
      graph_t *g = new graph_t;

      vector<tuple<shared_ptr<Vertex>, shared_ptr<Vertex>, int, int>> at_list; // (v, w, t, col_idx)
      int col_idx = 0;
      for (auto [col, val] : mas_result->solution)
      {
        int time = 0;
        int v_idx = 0;
        for (auto v : col->path)
        {
          if (v_idx >= 1 && v_idx <= col->path.size() - 2)
            at_list.push_back(make_tuple(v, col->path[v_idx + 1], time, col_idx));
          if (v->type == Vertex::Type::NORMAL)
            time += 1;
          v_idx += 1;
        }
        col_idx += 1;
      }
      g->n = at_list.size();
      g->edges = static_cast<set_t *>(calloc(g->n, sizeof(set_t)));
      for (int i = 0; i < g->n; i++)
        g->edges[i] = set_new(g->n);
      g->weights = static_cast<int *>(calloc(g->n, sizeof(int)));
      for (int i = 0; i < g->n; i++)
        g->weights[i] = (int)std::floor(cl_threshold * std::get<1>(mas_result->solution[std::get<3>(at_list[i])]));
      int at_idx1 = -1;
      for (auto [v1, w1, t1, col_idx1] : at_list)
      {
        at_idx1 += 1;
        int at_idx2 = -1;
        for (auto [v2, w2, t2, col_idx2] : at_list)
        {
          at_idx2 += 1;
          if (col_idx2 <= col_idx1)
          {
            continue;
          }
          bool is_incompatible = false;
          if (std::get<0>(mas_result->solution[col_idx1])->agent_id == std::get<0>(mas_result->solution[col_idx2])->agent_id)
          {
            if (t1 == t2)
              is_incompatible = true;
            else if (t1 <= t2 && t2 <= t1 + instance->param.time_spacing && v1 == v2 && v1 != w1)
              is_incompatible = true;
            else if (t1 <= t2 + 1 && t2 + 1 <= t1 + instance->param.time_spacing && v1 == w2 && v1 != w1)
              is_incompatible = true;
            else if (t2 <= t1 && t1 <= t2 + instance->param.time_spacing && v1 == v2 && v2 != w2)
              is_incompatible = true;
            else if (t2 <= t1 + 1 && t1 + 1 <= t2 + instance->param.time_spacing && v2 == w1 && v2 != w2)
              is_incompatible = true;
          }
          else
          {
            if (v1 == v2 && t1 == t2)
              is_incompatible = true;
            else if (w1 == w2 && t1 == t2)
              is_incompatible = true;
            else if (v1 == w2 && v2 == w1 && t1 == t2)
              is_incompatible = true;
            else if (std::abs(t1 - t2) <= instance->param.time_spacing && v1 == v2)
              is_incompatible = true;
            else if (std::abs(t1 - t2 - 1) <= instance->param.time_spacing && v1 == w2)
              is_incompatible = true;
            else if (std::abs(t1 + 1 - t2) <= instance->param.time_spacing && w1 == v2)
              is_incompatible = true;
            else if (std::abs(t1 - t2) <= instance->param.time_spacing && w1 == w2)
              is_incompatible = true;
          }
          if (is_incompatible)
            GRAPH_ADD_EDGE(g, at_idx1, at_idx2);
        }
      }
      // 1. separate violating cliques

      // Find all violating cliques
      int num_violated_cliques = clique_find_all(g, cl_threshold, 0, TRUE, opts);
      // println("num {} {}", num_violated_cliques, opts->clique_list_length);
      for (int i = 0; i < std::min(num_violated_cliques, opts->clique_list_length); i++)
      {
        auto clique = opts->clique_list[i];
        // println("val {}", graph_subgraph_weight(g, clique));

        Clique cl;
        int idx = -1;
        while ((idx = set_return_next(clique, idx)) >= 0)
        {
          auto [v, w, t, col_idx] = at_list[idx];
          auto [col, _] = mas_result->solution[col_idx];
          // println("col_idx {}  val {}  i {}  v {}  w {}  t {}", col_idx, g->weights[idx], col->agent_id, v->name(), w->name(), t);
          cl.iat_list.push_back(make_tuple(col->agent_id, v, w, t, col));
          // println("idx {}  val {}", idx, g->weights[idx])
          // auto [col, _] = mas_result->solution[idx];
          col->clique_conflict_coeff[cliques->size()] = 1;
          col->is_altered = true;
        }
        cut_added = true;
        num_new_cut += 1;
        num_cc_ctrs += 1;
        cliques->push_back(cl);
      }
      mas_solver->cliques = cliques;
      graph_free(g);
      free(opts->clique_list);
    }
    total_time_spent_on_sepa += get_elapsed(sepa_start);
    // Add conflict 4 sepa

    mas_solver->add_separated_cuts(local_pool);

#ifdef SHOW_COLGEN_RESULT
    println("      Cut-{}: obj {:.3f}  #row {}  #col {}  nnc {}  v {}  a {}  t {}  w {}  r {}  o {}  c {}  1 {}  2 {}  3 {}  4 {}  time {:.3f}", num_cut_iter, mas_result->obj_val, mas_result->num_row, mas_result->num_col, num_new_cut, num_vc_ctrs, num_ac_ctrs, num_tc_ctrs, num_wc_ctrs, num_rc_ctrs, num_oc_ctrs, num_cc_ctrs, num_1c_ctrs, num_2c_ctrs, num_3c_ctrs, num_4c_ctrs, get_elapsed(node_start));
    ; // 6 initial spaces
#endif

    num_cut_iter += 1;
  }

  // if (depth == 0) // && instance->param.algorithm == Param::Alg::bpc)
  // {
  //   println("reduced_cost_fixing");
  //   auto rcf_start = get_now();
  //   int num_fixed = 0;
  //   for (int i = 0; i < instance->param.num_agents; i++)
  //   {
  //     // i = 6;
  //     num_fixed += reduced_cost_fixing(instance->graph, i, mas_result, erased_arcs[i], ub);
  //     // break;
  //   }
  //   println("reduced_cost_fixing done in {:.3f}  #fixed {}", get_elapsed(rcf_start), num_fixed);
  // }

  total_time_spent = get_elapsed(node_start);

  num_explored++;
}

void BPCNode::print_result(int num_unexplored, bool node_feasible, bool node_suboptimal, bool node_fractional, double global_ub)
{
  if (parent_lb == DOUBLE_NULL)
    parent_lb = -TRIVIAL_BOUND;
  if (ub == DOUBLE_NULL)
    ub = TRIVIAL_BOUND;
  if (global_ub == DOUBLE_NULL)
    global_ub = TRIVIAL_BOUND;
  print("    Node-{}: id {}  depth {}  #unexplored  {}  #total_cols {}  #local_cols {}  parent_lb {:.3f}  node_lb {:.3f}  global_ub {:.3f}  #iter {}  avg_sub_time {:.3f}  sub_time {:.3f}  mas_time {:.3f}  time {:.3f}", num_explored, id, depth, num_unexplored, pool->size(), local_pool.size(), parent_lb, lb, global_ub, num_iter, avg_time_spent_on_sp, total_time_spent_on_sub, total_time_spent_on_mas, total_time_spent); // 4 initial spaces
  print("  status {}", node_feasible ? (node_suboptimal ? "suboptimal" : (node_fractional ? "fractional" : "integer")) : "infeasible");
  println("");
}

bool BPCNode::is_feasible()
{
  return std::none_of(base_columns.begin(), base_columns.end(),
                      [](const auto &cc)
                      { return std::get<0>(cc)->is_dummy; });
}

bool BPCNode::has_fractional_solution()
{
  bool fractional_found = false;
  for (auto &[col, value] : base_columns)
  {
    double fractional_part = value - ((int)value);
    if (fractional_part < 0)
      throw;
    if (fractional_part > INTEGRALITY_EPS && fractional_part < 1 - INTEGRALITY_EPS)
    {
      fractional_found = true;
      break;
    }
  }
  return fractional_found;
}


MAPFSol BPCNode::primal_heuristic()
{
  double obj_val = 0;
  vector<vector<Coord>> paths;

  sy_map<tuple<shared_ptr<Vertex>, int>, bool> vt_occupied;

  for (auto pair : agent_order)
  {
    auto agent_id = std::get<0>(pair);

    auto [v_path, _, cost] = solve_spp_teg_astar(instance->graph,
                                                 agent_id,
                                                 instance->param.num_agents,
                                                 instance->param.time_spacing,
                                                 mas_result,
                                                 erased_arcs[agent_id],
                                                 vt_occupied);

    if (cost == DUMMY_COLUMN_COST)
    {
      obj_val = DUMMY_COLUMN_COST;
      paths = vector<vector<Coord>>();
      return MAPFSol(instance, obj_val, paths);
    }
    vector<Coord> path;
    int time = -1;
    for (auto v : v_path)
    {
      if (v->type == Vertex::Type::NORMAL)
      {
        time += 1;
        path.push_back(v->coord);
        for (int tt = -instance->param.time_spacing; tt <= instance->param.time_spacing; tt++)
          vt_occupied[make_tuple(v, time + tt)] = true;
      }
    }
    paths.push_back(path);
    obj_val += cost;
  }
  return MAPFSol(instance, obj_val, paths);
}
