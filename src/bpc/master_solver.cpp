#include "master_solver.h"
#include "sy_log.h"


MasterSolver::~MasterSolver()
{
  env.end();
}

void MasterSolver::reset_model()
{
  num_added_columns = 0;
  is_just_reset = true;
  std::stringstream temp_name;

  env = IloEnv();
  model = IloModel(env);
  vars = IloNumVarArray(env);
  pat_vars = IloNumVarArray(env);

  vertex_conflict_ctrs = sy_map<tuple<int, int>, IloRange>();
  vertex_conflict_is_needed = sy_map<tuple<int, int>, bool>();
  vertex_conflict_is_added = sy_map<tuple<int, int>, bool>();
  arc_conflict_ctrs = sy_map<tuple<int, int>, IloRange>();
  arc_conflict_is_needed = sy_map<tuple<int, int>, bool>();
  arc_conflict_is_added = sy_map<tuple<int, int>, bool>();
  time_conflict_ctrs = sy_map<tuple<int, int, int, int>, IloRange>();
  time_conflict_is_needed = sy_map<tuple<int, int, int, int>, bool>();
  time_conflict_is_added = sy_map<tuple<int, int, int, int>, bool>();
  wait_conflict_ctrs = sy_map<tuple<int, int>, IloRange>();
  wait_conflict_is_needed = sy_map<tuple<int, int>, bool>();
  wait_conflict_is_added = sy_map<tuple<int, int>, bool>();
  range_conflict_ctrs = sy_map<tuple<int, int>, IloRange>();
  range_conflict_is_needed = sy_map<tuple<int, int>, bool>();
  range_conflict_is_added = sy_map<tuple<int, int>, bool>();
  old_conflict_ctrs = sy_map<tuple<int, int, int>, IloRange>();
  old_conflict_is_needed = sy_map<tuple<int, int, int>, bool>();
  old_conflict_is_added = sy_map<tuple<int, int, int>, bool>();
  clique_conflict_ctrs = sy_map<int, IloRange>(); // cc
  clique_conflict_is_needed = sy_map<int, bool>();
  clique_conflict_is_added = sy_map<int, bool>();
  tvc1_conflict_ctrs = sy_map<tuple<int, int, int, int>, IloRange>();
  tvc1_conflict_is_needed = sy_map<tuple<int, int, int, int>, bool>();
  tvc1_conflict_is_added = sy_map<tuple<int, int, int, int>, bool>();
  tvc2_conflict_ctrs = sy_map<tuple<int, int, int, int>, IloRange>();
  tvc2_conflict_is_needed = sy_map<tuple<int, int, int, int>, bool>();
  tvc2_conflict_is_added = sy_map<tuple<int, int, int, int>, bool>();
  tvc3_conflict_ctrs = sy_map<tuple<int, int, int, int>, IloRange>();
  tvc3_conflict_is_needed = sy_map<tuple<int, int, int, int>, bool>();
  tvc3_conflict_is_added = sy_map<tuple<int, int, int, int>, bool>();
  tvc4_conflict_ctrs = sy_map<tuple<int, int, int, int>, IloRange>();
  tvc4_conflict_is_needed = sy_map<tuple<int, int, int, int>, bool>();
  tvc4_conflict_is_added = sy_map<tuple<int, int, int, int>, bool>();
  /// Add conflict 5 clear

  /* Declare constraints */
  obj = IloMinimize(env, 0.0); // columns will be added later

  /* Generate constraints */
  // Cardinality constraints
  for (int agent_id = 0; agent_id < instance->param.num_agents; agent_id++)
  {
    temp_name.str("");
    temp_name << "cardi_" << agent_id;
    cardi_ctrs[agent_id] = IloRange(env, 1.0, 1.0, temp_name.str().c_str());
  }
}

void MasterSolver::generate_model(const ColumnPool &pool, bool solve_linear, string *model_dump)
{
  std::stringstream temp_name;

  // Pattern variables
  int num_patterns = 0;
  int idx = -1;
  for (auto col : pool)
  {
    idx += 1;
    if (num_patterns < num_added_columns)
    {
      num_patterns += 1;
      continue;
    }
    else // for new columns (num_patterns >= num_added_columns)
    {
      // extend time conflict
      if (!col->is_dummy)
      {
        for (auto &[key, _] : time_conflict_ctrs)
        {
          auto &[agent_id, vertex_id, time, delt] = key;
          if (col->agent_id == agent_id && time + 1 <= col->path.size() - 1 && col->path[time + 1]->id == vertex_id)
            col->time_conflict_coeff[key] = 1;

          else if (col->agent_id != agent_id)
          {
            if (time + delt + 1 <= col->path.size() - 1)
            {
              if (col->path[time + delt + 1]->id == vertex_id)
                col->time_conflict_coeff[key] = 1;
            }
          }
        }
        sy_map<tuple<int, int, int>, double> key_coeff;
        for (auto &[key, _] : old_conflict_ctrs)
        {
          auto &[agent_id, vertex_id, time] = key;
          if (col->agent_id == agent_id && time + 1 <= col->path.size() - 1 && col->path[time + 1]->id == vertex_id)
            key_coeff[key] += 1;
          else if (col->agent_id != agent_id)
          {
            for (int delt = 0; delt <= instance->param.time_spacing; delt++)
              if (time + delt + 1 <= col->path.size() - 1 && col->path[time + delt + 1]->id == vertex_id)
                key_coeff[key] = 1 / ((double)1 + instance->param.time_spacing);
          }
          col->old_conflict_coeff[key] = key_coeff[key];
        }
      }
    }

    // Objective function coefficient
    IloNumColumn ilo_col = obj(col->obj_coeff);

    // Cardinality constraints coefficient
    for (int agent_id = 0; agent_id < instance->param.num_agents; agent_id++)
      ilo_col += cardi_ctrs[agent_id](col->cardi_coeff[agent_id]);

    // Vertex conflict constraints coefficient
    if (instance->param.ctrs_string.find('V') != string::npos || instance->param.formul_string.find('V') != string::npos)
      for (auto [key, coeff] : col->vertex_conflict_coeff)
      {
        if (!vertex_conflict_ctrs.contains(key))
        {
          temp_name << fmt::format("vc_{}_{}", instance->graph.get_vertex(std::get<0>(key))->name(), std::get<1>(key));
          vertex_conflict_ctrs[key] = IloRange(env, 0, 1.0, temp_name.str().c_str());
          temp_name.str("");
        }
        else
          vertex_conflict_is_needed[key] = true;
        
        ilo_col += vertex_conflict_ctrs[key](coeff);
      }
    // Arc conflict constraints coefficient
    if (instance->param.ctrs_string.find('A') != string::npos || instance->param.formul_string.find('A') != string::npos)
      for (auto [key, coeff] : col->arc_conflict_coeff)
      {
        if (!arc_conflict_ctrs.contains(key))
        {
          temp_name << fmt::format("ac_{}_{}", instance->graph.get_arc(std::get<0>(key))->name(), std::get<1>(key));
          arc_conflict_ctrs[key] = IloRange(env, 0, 1.0, temp_name.str().c_str());
          temp_name.str("");
        }
        else
          arc_conflict_is_needed[key] = true;
        ilo_col += arc_conflict_ctrs[key](coeff);
      }
    // Time conflict constraints coefficient
    if (instance->param.ctrs_string.find('T') != string::npos)
      for (auto [key, coeff] : col->time_conflict_coeff)
      {
        auto [a, v, t, d] = key;
        if (!time_conflict_ctrs.contains(key))
        {
          temp_name << fmt::format("tc_{}_{}_{}_{}", a, instance->graph.get_vertex(v)->name(), t, d);
          time_conflict_ctrs[key] = IloRange(env, 0, 1.0, temp_name.str().c_str());
          temp_name.str("");
        }
        else
          time_conflict_is_needed[key] = true;
        ilo_col += time_conflict_ctrs[key](coeff);
      }
    // Wait conflict constraints coefficient
    if (instance->param.ctrs_string.find('W') != string::npos || instance->param.formul_string.find('W') != string::npos)
      for (auto [key, coeff] : col->wait_conflict_coeff)
      {
        if (!wait_conflict_ctrs.contains(key))
        {
          temp_name << fmt::format("wc_{}_{}", instance->graph.get_vertex(std::get<0>(key))->name(), std::get<1>(key));
          wait_conflict_ctrs[key] = IloRange(env, 0, 1.0, temp_name.str().c_str());
          temp_name.str("");
        }
        else
          wait_conflict_is_needed[key] = true;
        ilo_col += wait_conflict_ctrs[key](coeff);
      }
    // Range conflict constraints coefficient
    if (instance->param.ctrs_string.find('R') != string::npos || instance->param.formul_string.find('R') != string::npos)
      for (auto [key, coeff] : col->range_conflict_coeff)
      {
        if (!range_conflict_ctrs.contains(key))
        {
          temp_name << fmt::format("rc_{}_{}", instance->graph.get_vertex(std::get<0>(key))->name(), std::get<1>(key));
          range_conflict_ctrs[key] = IloRange(env, 0, 1.0, temp_name.str().c_str());
          temp_name.str("");
        }
        else
          range_conflict_is_needed[key] = true;
        ilo_col += range_conflict_ctrs[key](coeff);
      }
    // old conflict constraints coefficient
    if (instance->param.ctrs_string.find('O') != string::npos)
      for (auto [key, coeff] : col->old_conflict_coeff)
      {
        auto [a, v, t] = key;
        if (!old_conflict_ctrs.contains(key))
        {
          temp_name << fmt::format("oc_{}_{}_{}", a, instance->graph.get_vertex(v)->name(), t);
          old_conflict_ctrs[key] = IloRange(env, 0, 1.0, temp_name.str().c_str());
          temp_name.str("");
        }
        else
          old_conflict_is_needed[key] = true;
        ilo_col += old_conflict_ctrs[key](coeff);
      }
    // clique conflict constraints coefficient
    if (instance->param.ctrs_string.find('C') != string::npos)
      for (auto [key, coeff] : col->clique_conflict_coeff)
      {
        auto idx = key;
        if (!clique_conflict_ctrs.contains(key))
        {
          temp_name << fmt::format("cc_{}", idx);
          clique_conflict_ctrs[key] = IloRange(env, 0, 1.0, temp_name.str().c_str());
          temp_name.str("");
        }
        else
          clique_conflict_is_needed[key] = true;
        ilo_col += clique_conflict_ctrs[key](coeff);
      }
    // tvc1 conflict constraints coefficient
    if (instance->param.ctrs_string.find('1') != string::npos)
      for (auto [key, coeff] : col->tvc1_conflict_coeff)
      {
        auto [i, v_id, w_id, t] = key;
        if (!tvc1_conflict_ctrs.contains(key))
        {
          temp_name << fmt::format("1c_{}_{}_{}_{}", i, v_id, w_id, t);
          tvc1_conflict_ctrs[key] = IloRange(env, 0, 1.0, temp_name.str().c_str());
          temp_name.str("");
        }
        else
          tvc1_conflict_is_needed[key] = true;
        ilo_col += tvc1_conflict_ctrs[key](coeff);
      }
    // tvc2 conflict constraints coefficient
    if (instance->param.ctrs_string.find('2') != string::npos)
      for (auto [key, coeff] : col->tvc2_conflict_coeff)
      {
        auto [i, v_id, w_id, t] = key;
        if (!tvc2_conflict_ctrs.contains(key))
        {
          temp_name << fmt::format("2c_{}_{}_{}_{}", i, v_id, w_id, t);
          tvc2_conflict_ctrs[key] = IloRange(env, 0, 1.0, temp_name.str().c_str());
          temp_name.str("");
        }
        else
          tvc2_conflict_is_needed[key] = true;
        ilo_col += tvc2_conflict_ctrs[key](coeff);
      }
    // tvc3 conflict constraints coefficient
    if (instance->param.ctrs_string.find('3') != string::npos)
      for (auto [key, coeff] : col->tvc3_conflict_coeff)
      {
        auto [i, v_id, w_id, t] = key;
        if (!tvc3_conflict_ctrs.contains(key))
        {
          temp_name << fmt::format("3c_{}_{}_{}_{}", i, v_id, w_id, t);
          tvc3_conflict_ctrs[key] = IloRange(env, 0, 1.0, temp_name.str().c_str());
          temp_name.str("");
          // tvc3_conflict_is_needed[key] = true;
        }
        else
          tvc3_conflict_is_needed[key] = true;
        ilo_col += tvc3_conflict_ctrs[key](coeff);
      }
    // tvc4 conflict constraints coefficient
    if (instance->param.ctrs_string.find('4') != string::npos)
      for (auto [key, coeff] : col->tvc4_conflict_coeff)
      {
        auto [i, v_id, w_id, t] = key;
        if (!tvc4_conflict_ctrs.contains(key))
        {
          temp_name << fmt::format("4c_{}_{}_{}_{}", i, v_id, w_id, t);
          tvc4_conflict_ctrs[key] = IloRange(env, 0, 1.0, temp_name.str().c_str());
          temp_name.str("");
        }
        else
          tvc4_conflict_is_needed[key] = true;
        ilo_col += tvc4_conflict_ctrs[key](coeff);
      }
    // Add conflict 6 gene

    // Generate variable (column-wise)
    temp_name << (col->is_dummy ? "dummy" : "p" + to_string(num_patterns));
    IloNumVar var(ilo_col, 0, IloInfinity, (solve_linear ? IloNumVar::Float : IloNumVar::Int), temp_name.str().c_str());
    temp_name.str("");

    pat_vars.add(var);
    num_added_columns += 1;
    col->is_altered = false;
    num_patterns += 1;
  }

  /* Generate model and solve */
  // objective
  if (is_just_reset)
    model.add(obj);

  // constraints
  auto ctrs = IloRangeArray(env);
  if (is_just_reset)
  {
    for (int agent_id = 0; agent_id < instance->param.num_agents; agent_id++)
      ctrs.add(cardi_ctrs[agent_id]);
    is_just_reset = false;
  }

  for (auto &[key, ctr] : vertex_conflict_ctrs)
    if (vertex_conflict_is_needed.contains(key) && !vertex_conflict_is_added.contains(key))
    {
      ctrs.add(ctr);
      vertex_conflict_is_added[key] = true;
    }
  for (auto &[key, ctr] : arc_conflict_ctrs)
    if (arc_conflict_is_needed.contains(key) && !arc_conflict_is_added.contains(key))
    {
      ctrs.add(ctr);
      arc_conflict_is_added[key] = true;
    }
  for (auto &[key, ctr] : time_conflict_ctrs)
    if (time_conflict_is_needed.contains(key) && !time_conflict_is_added.contains(key))
    {
      ctrs.add(ctr);
      time_conflict_is_added[key] = true;
    }
  for (auto &[key, ctr] : wait_conflict_ctrs)
    if (wait_conflict_is_needed.contains(key) && !wait_conflict_is_added.contains(key))
    {
      ctrs.add(ctr);
      wait_conflict_is_added[key] = true;
    }
  for (auto &[key, ctr] : range_conflict_ctrs)
    if (range_conflict_is_needed.contains(key) && !range_conflict_is_added.contains(key))
    {
      ctrs.add(ctr);
      range_conflict_is_added[key] = true;
    }
  for (auto &[key, ctr] : old_conflict_ctrs)
    if (old_conflict_is_needed.contains(key) && !old_conflict_is_added.contains(key))
    {
      ctrs.add(ctr);
      old_conflict_is_added[key] = true;
    }
  for (auto &[key, ctr] : clique_conflict_ctrs)
    if (clique_conflict_is_needed.contains(key) && !clique_conflict_is_added.contains(key))
    {
      ctrs.add(ctr);
      clique_conflict_is_added[key] = true;
    }
  for (auto &[key, ctr] : tvc1_conflict_ctrs)
    if (tvc1_conflict_is_needed.contains(key) && !tvc1_conflict_is_added.contains(key))
    {
      ctrs.add(ctr);
      tvc1_conflict_is_added[key] = true;
    }
  for (auto &[key, ctr] : tvc2_conflict_ctrs)
    if (tvc2_conflict_is_needed.contains(key) && !tvc2_conflict_is_added.contains(key))
    {
      ctrs.add(ctr);
      tvc2_conflict_is_added[key] = true;
    }
  for (auto &[key, ctr] : tvc3_conflict_ctrs)
    if (tvc3_conflict_is_needed.contains(key) && !tvc3_conflict_is_added.contains(key))
    {
      ctrs.add(ctr);
      tvc3_conflict_is_added[key] = true;
    }
  for (auto &[key, ctr] : tvc4_conflict_ctrs)
    if (tvc4_conflict_is_needed.contains(key) && !tvc4_conflict_is_added.contains(key))
    {
      ctrs.add(ctr);
      tvc4_conflict_is_added[key] = true;
    }
  model.add(ctrs);
  // Add conflict 7 add coef

  cplex = IloCplex(model);

  if (model_dump)
  {
    string output_name = *model_dump + ".lp";
    try
    {
      println("model dump at {}", output_name);
      cplex.exportModel(output_name.c_str());
    }
    catch (IloException &e)
    {
      std::cerr << "Export IloException: " << e << std::endl;
    }
  }
}

void MasterSolver::add_separated_cuts(const ColumnPool &pool)
{
  if (instance->param.ctrs_string.find('1') != string::npos)
  {
    for (auto [key, needed]: tvc1_conflict_is_needed)
    {
      if (!needed)
        continue;

      IloNumVarArray vars(env);
      IloNumArray vals(env);
      int idx = -1;
      for (auto col : pool) {
        idx += 1;
        if (!col->tvc1_conflict_coeff.contains(key))
          continue;
        vars.add(pat_vars[idx]);
        vals.add(col->tvc1_conflict_coeff[key]);
      }
      IloRange ilocut(env, 0, 1.0);
      ilocut.setLinearCoefs(vars, vals);
      model.add(ilocut);
      tvc1_conflict_ctrs[key] = ilocut;
      vars.end();
      vals.end();
      tvc1_conflict_is_added[key] = true;
    }
    tvc1_conflict_is_needed = sy_map<tuple<int, int, int, int>, bool>();
  }
  if (instance->param.ctrs_string.find('2') != string::npos)
  {
    for (auto [key, needed]: tvc2_conflict_is_needed)
    {
      if (!needed)
        continue;

      IloNumVarArray vars(env);
      IloNumArray vals(env);
      int idx = -1;
      for (auto col : pool) {
        idx += 1;
        if (!col->tvc2_conflict_coeff.contains(key))
          continue;
        vars.add(pat_vars[idx]);
        vals.add(col->tvc2_conflict_coeff[key]);
      }
      IloRange ilocut(env, 0, 1.0);
      ilocut.setLinearCoefs(vars, vals);
      model.add(ilocut);
      tvc2_conflict_ctrs[key] = ilocut;
      vars.end();
      vals.end();
      tvc2_conflict_is_added[key] = true;
    }
    tvc2_conflict_is_needed = sy_map<tuple<int, int, int, int>, bool>();
  }
  if (instance->param.ctrs_string.find('3') != string::npos)
  {
    for (auto [key, needed]: tvc3_conflict_is_needed)
    {
      if (!needed)
        continue;

      IloNumVarArray vars(env);
      IloNumArray vals(env);
      int idx = -1;
      for (auto col : pool) {
        idx += 1;
        if (!col->tvc3_conflict_coeff.contains(key))
          continue;
        vars.add(pat_vars[idx]);
        vals.add(col->tvc3_conflict_coeff[key]);
      }
      IloRange ilocut(env, 0, 1.0);
      ilocut.setLinearCoefs(vars, vals);
      model.add(ilocut);
      tvc3_conflict_ctrs[key] = ilocut;
      vars.end();
      vals.end();
      tvc3_conflict_is_added[key] = true;
    }
    tvc3_conflict_is_needed = sy_map<tuple<int, int, int, int>, bool>();
  }
  if (instance->param.ctrs_string.find('4') != string::npos)
  {
    for (auto [key, needed]: tvc4_conflict_is_needed)
    {
      if (!needed)
        continue;

      IloNumVarArray vars(env);
      IloNumArray vals(env);
      int idx = -1;
      for (auto col : pool) {
        idx += 1;
        if (!col->tvc4_conflict_coeff.contains(key))
          continue;
        vars.add(pat_vars[idx]);
        vals.add(col->tvc4_conflict_coeff[key]);
      }
      IloRange ilocut(env, 0, 1.0);
      ilocut.setLinearCoefs(vars, vals);
      model.add(ilocut);
      tvc4_conflict_ctrs[key] = ilocut;
      vars.end();
      vals.end();
      tvc4_conflict_is_added[key] = true;
    }
    tvc4_conflict_is_needed = sy_map<tuple<int, int, int, int>, bool>();
  }
}

void MasterSolver::solve(const ColumnPool &pool, bool solve_linear, string *model_dump)
{
  auto mas_solve_start = get_now();
  generate_model(pool, solve_linear, model_dump);
#ifndef SHOW_CPLEX_LOG
  cplex.setOut(env.getNullStream());
#endif

  auto solved = false;
  try
  {
    solved = cplex.solve();
  }
  catch (IloException &e)
  {
    std::cerr << "Solve IloException: " << e << std::endl;
    throw;
  }

  if (!solved)
  {
    println(">> Integer feasible solution not found (bpc_master)");
  }
}

/*
Wrapper method for MasterSolver::solve(pool, solve_linear)
Summarize LP results - objval, duals, solution - and return
*/
shared_ptr<MasterResult> MasterSolver::solve_lp(const ColumnPool &pool, string *model_dump)
{
#ifdef SHOW_MASSUB_RESULT
  auto mas_lp_start = get_now(); 
#endif
  solve(pool, true, model_dump);
  double obj_val = cplex.getObjValue();

#ifdef SHOW_MASSUB_RESULT
  println("            MP: #cols {}  #rows {}  stat {}  time {}", cplex.getNcols(), cplex.getNrows(), cplex.getStatus(), get_elapsed(mas_lp_start)); // 12 initial spaces
#endif

  sy_map<int, sy_map<tuple<shared_ptr<Arc>, int>, double>> i_at_cost_altered;
  sy_map<tuple<shared_ptr<Arc>, int>, double> at_cost_altered;

  /* Get duals of constraints */
  auto cardi_duals = i_values();
  for (int agent_id = 0; agent_id < instance->param.num_agents; agent_id++)
  {
    double value = cplex.getDual(cardi_ctrs[agent_id]);
    if (-DUAL_VALUE_EPS <= value && value <= DUAL_VALUE_EPS)
      continue;
    cardi_duals[agent_id] = value;
  }
  // Vertex conflict duals
  auto vertex_conflict_duals = vt_values();
  if (instance->param.ctrs_string.find('V') != std::string::npos || instance->param.formul_string.find('V') != string::npos)
  {
    for (auto &[key, _] : vertex_conflict_is_needed)
    {
      double value = cplex.getDual(vertex_conflict_ctrs[key]);
      if (value < 0)
      {
        vertex_conflict_duals[key] = value;

        auto &[v_id, t] = key;
        auto v = instance->graph.get_vertex(v_id);
        for (auto a : instance->graph.get_out_arcs(v))
        {
          at_cost_altered[make_tuple(a, t)] -= value;
        }
      }
    }
  }
  // Arc conflict duals
  auto arc_conflict_duals = at_values();
  if (instance->param.ctrs_string.find('A') != std::string::npos || instance->param.formul_string.find('A') != string::npos)
  {
    for (auto &[key, _] : arc_conflict_is_needed)
    {
      double value = cplex.getDual(arc_conflict_ctrs[key]);

      if (value < 0)
      {
        arc_conflict_duals[key] = value;

        auto &[a_id, t] = key;
        auto a = instance->graph.get_arc(a_id);
        auto ar = instance->graph.get_arc(a->to, a->from);
        at_cost_altered[make_tuple(a, t)] -= value;
        at_cost_altered[make_tuple(ar, t)] -= value;
      }
    }
  }
  // Time conflict duals
  auto time_conflict_duals = ivtt_values();
  if (instance->param.ctrs_string.find('T') != std::string::npos)
  {
    for (auto &[key, _] : time_conflict_ctrs)
    {
      double value = cplex.getDual(time_conflict_ctrs[key]);
      time_conflict_duals[key] = value; 
      if (value < 0)
      {
        auto &[agent_id, v_id, t, delt] = key;
        auto v = instance->graph.get_vertex(v_id);
        for (auto a : instance->graph.get_out_arcs(v))
        {
          i_at_cost_altered[agent_id][make_tuple(a, t)] -= value;
          for (int j = 0; j < instance->param.num_agents; j++)
            if (j != agent_id)
              i_at_cost_altered[j][make_tuple(a, t + delt)] -= value;
        }
      }
    }
  }
  // Wait conflict duals
  auto wait_conflict_duals = vt_values();
  if (instance->param.ctrs_string.find('W') != std::string::npos || instance->param.formul_string.find('W') != string::npos)
  {
    for (auto &[key, _] : wait_conflict_is_needed)
    {
      double value = cplex.getDual(wait_conflict_ctrs[key]);

      if (value < 0)
      {
        wait_conflict_duals[key] = value;

        auto &[v_id, t] = key;
        auto v = instance->graph.get_vertex(v_id);
        for (auto a : instance->graph.get_out_arcs(v))
        {
          at_cost_altered[make_tuple(a, t)] -= value;
        }
        for (auto a : instance->graph.get_inin_arcs(v))
        {
          at_cost_altered[make_tuple(a, t)] -= value;
        }
      }
    }
  }
  // Range conflict duals
  auto range_conflict_duals = vt_values();
  if (instance->param.ctrs_string.find('R') != std::string::npos || instance->param.formul_string.find('R') != string::npos)
  {
    for (auto &[key, _] : range_conflict_is_needed)
    {
      double value = cplex.getDual(range_conflict_ctrs[key]);

      if (value < 0)
      {
        range_conflict_duals[key] = value;

        auto &[v_id, t] = key;
        auto v = instance->graph.get_vertex(v_id);
        for (auto a : instance->graph.get_out_arcs(v))
        {
          if (a->from != a->to)
            for (int delt = 0; delt <= instance->param.time_spacing; delt++)
              at_cost_altered[make_tuple(a, t + delt)] -= value;
          else
            at_cost_altered[make_tuple(a, t + instance->param.time_spacing)] -= value;
        }
      }
    }
  }
  // Old conflict duals
  auto old_conflict_duals = ivt_values();
  if (instance->param.ctrs_string.find('O') != std::string::npos)
  {
    for (auto &[key, _] : old_conflict_ctrs)
    {
      double value = cplex.getDual(old_conflict_ctrs[key]);
      old_conflict_duals[key] = value;
      if (value < 0)
      {
        auto &[agent_id, v_id, t] = key;
        auto v = instance->graph.get_vertex(v_id);
        for (auto a : instance->graph.get_out_arcs(v))
        {
          i_at_cost_altered[agent_id][make_tuple(a, t)] -= value;
          for (int j = 0; j < instance->param.num_agents; j++)
          {
            if (j != agent_id)
            {
              for (int delt = 0; delt <= instance->param.time_spacing; delt++)
              {
                i_at_cost_altered[j][make_tuple(a, t + delt)] -= (double)value / (1 + instance->param.time_spacing);
              }
            }
          }
        }
      }
    }
  }
  // clique conflict duals
  auto clique_conflict_duals = sy_map<int, double>();
  if (instance->param.ctrs_string.find('C') != std::string::npos)
  {
    for (auto &[key, _] : clique_conflict_ctrs)
    {
      double value = cplex.getDual(clique_conflict_ctrs[key]);
      if (value < 0)
      {
        clique_conflict_duals[key] = value;
        auto clique = cliques->at(key);
        for (auto [agent_id, v, w, t, col] : clique.iat_list)
        {
          auto a = instance->graph.get_arc(v, w);
          i_at_cost_altered[agent_id][make_tuple(a, t)] -= value;
        }
      }
    }
  }
  // Two-vertex-clique-1 conflict duals
  auto tvc1_conflict_duals = ivwt_values();
  if (instance->param.ctrs_string.find('1') != std::string::npos)
  {
    for (auto &[key, _] : tvc1_conflict_is_added)
    {
      double value = cplex.getDual(tvc1_conflict_ctrs[key]);
      tvc1_conflict_duals[key] = value;
      if (value < 0)
      {
        auto &[agent_id, v_id, w_id, t] = key;
        if (v_id <= w_id)
          exit(1);
        auto v = instance->graph.get_vertex(v_id);
        auto w = instance->graph.get_vertex(w_id);
        for (auto a : instance->graph.get_out_arcs(v))
        {
          i_at_cost_altered[agent_id][make_tuple(a, t)] -= value;
        }
        for (auto a : instance->graph.get_out_arcs(w))
        {
          i_at_cost_altered[agent_id][make_tuple(a, t)] -= value;
        }
        for (auto a : instance->graph.get_inin_arcs(v))
        {
          if (a->from != v && a->from != w)
            i_at_cost_altered[agent_id][make_tuple(a, t)] -= value;
        }
        for (auto a : instance->graph.get_inin_arcs(w))
        {
          if (a->from != v && a->from != w)
            i_at_cost_altered[agent_id][make_tuple(a, t)] -= value;
        }

        for (int j = 0; j < instance->param.num_agents; j++)
        {
          if (j != agent_id)
          {
            for (int delt = 0; delt <= instance->param.time_spacing - 1; delt++)
            {
              i_at_cost_altered[j][make_tuple(instance->graph.get_arc(v, w), t + delt)] -= value;
              i_at_cost_altered[j][make_tuple(instance->graph.get_arc(w, v), t + delt)] -= value;
            }
          }
        }
      }
    }
  }
  // Two-vertex-clique-2 conflict duals
  auto tvc2_conflict_duals = ivwt_values();
  if (instance->param.ctrs_string.find('2') != std::string::npos)
  { 
    for (auto &[key, _] : tvc2_conflict_is_added)
    {
      double value = cplex.getDual(tvc2_conflict_ctrs[key]);
      tvc2_conflict_duals[key] = value;
      if (value < 0)
      {
        auto &[agent_id, v_id, w_id, t] = key;
        auto v = instance->graph.get_vertex(v_id);
        auto w = instance->graph.get_vertex(w_id);
        for (auto a : instance->graph.get_out_arcs(v))
        {
          i_at_cost_altered[agent_id][make_tuple(a, t)] -= value;
        }
        for (auto a : instance->graph.get_out_arcs(w))
        {
          i_at_cost_altered[agent_id][make_tuple(a, t)] -= value;
        }
        for (auto a : instance->graph.get_inin_arcs(v))
        {
          if (a->from != v && a->from != w)
            i_at_cost_altered[agent_id][make_tuple(a, t)] -= value;
        }
        for (auto a : instance->graph.get_inin_arcs(w))
        {
          if (a->from != v && a->from != w)
            i_at_cost_altered[agent_id][make_tuple(a, t)] -= value;
        }

        for (int j = 0; j < instance->param.num_agents; j++)
        {
          if (j != agent_id)
          {
            for (int delt = 0; delt <= instance->param.time_spacing - 1; delt++)
            {
              i_at_cost_altered[j][make_tuple(instance->graph.get_arc(v, w), t - delt)] -= value;
              i_at_cost_altered[j][make_tuple(instance->graph.get_arc(w, v), t - delt)] -= value;
            }
          }
        }
      }
    }
  }
  // Two-vertex-clique-3 conflict duals
  auto tvc3_conflict_duals = ivwt_values();
  if (instance->param.ctrs_string.find('3') != std::string::npos)
  {
    for (auto &[key, _] : tvc3_conflict_is_added)
    {
      double value = cplex.getDual(tvc3_conflict_ctrs[key]);
      tvc3_conflict_duals[key] = value;
      if (value < 0)
      {
        auto &[agent_id, v_id, w_id, t] = key;
        auto v = instance->graph.get_vertex(v_id);
        auto w = instance->graph.get_vertex(w_id);
        for (auto a : instance->graph.get_out_arcs(v))
        {
          i_at_cost_altered[agent_id][make_tuple(a, t)] -= value;
        }
        for (auto a : instance->graph.get_out_arcs(w))
        {
          if (a->to == v || a->to == w)
            i_at_cost_altered[agent_id][make_tuple(a, t)] -= value;
        }
        for (auto a : instance->graph.get_inin_arcs(v))
        {
          if (a->from != v && a->from != w)
            i_at_cost_altered[agent_id][make_tuple(a, t)] -= value;
        }
        for (auto a : instance->graph.get_inin_arcs(w))
        {
          if (a->from != v && a->from != w)
            i_at_cost_altered[agent_id][make_tuple(a, t)] -= value;
        }

        for (int j = 0; j < instance->param.num_agents; j++)
        {
          if (j != agent_id)
          {
            i_at_cost_altered[j][make_tuple(instance->graph.get_arc(v, w), t)] -= value;
            for (int delt = 1; delt <= instance->param.time_spacing - 1; delt++)
            {
              i_at_cost_altered[j][make_tuple(instance->graph.get_arc(v, w), t + delt)] -= value;
              i_at_cost_altered[j][make_tuple(instance->graph.get_arc(w, v), t + delt)] -= value;
            }
            i_at_cost_altered[j][make_tuple(instance->graph.get_arc(v, w), t + instance->param.time_spacing)] -= value;
          }
        }
      }
    }
  }
  // Two-vertex-clique-4 conflict duals
  auto tvc4_conflict_duals = ivwt_values();
  if (instance->param.ctrs_string.find('4') != std::string::npos)
  {
    for (auto &[key, _] : tvc4_conflict_is_added)
    {
      double value = cplex.getDual(tvc4_conflict_ctrs[key]);
      tvc4_conflict_duals[key] = value;
      if (value < 0)
      {
        auto &[agent_id, v_id, w_id, t] = key;
        auto v = instance->graph.get_vertex(v_id);
        auto w = instance->graph.get_vertex(w_id);
        for (auto a : instance->graph.get_out_arcs(v))
        {
          i_at_cost_altered[agent_id][make_tuple(a, t)] -= value;
        }
        for (auto a : instance->graph.get_out_arcs(w))
        {
          i_at_cost_altered[agent_id][make_tuple(a, t)] -= value;
        }
        for (auto a : instance->graph.get_inin_arcs(w))
        {
          if (a->from != v && a->from != w)
            i_at_cost_altered[agent_id][make_tuple(a, t)] -= value;
        }

        for (int j = 0; j < instance->param.num_agents; j++)
        {
          if (j != agent_id)
          {
            i_at_cost_altered[j][make_tuple(instance->graph.get_arc(v, w), t)] -= value;
            for (int delt = 1; delt <= instance->param.time_spacing - 1; delt++)
            {
              i_at_cost_altered[j][make_tuple(instance->graph.get_arc(v, w), t - delt)] -= value;
              i_at_cost_altered[j][make_tuple(instance->graph.get_arc(w, v), t - delt)] -= value;
            }
            i_at_cost_altered[j][make_tuple(instance->graph.get_arc(v, w), t - instance->param.time_spacing)] -= value;
          }
        }
      }
    }
  }
  // Add conflict 8 get dual

  IloNumArray values(env);
  cplex.getValues(values, pat_vars);

  auto solution = vector<tuple<shared_ptr<Column>, double>>();
  auto reserved = sy_map<tuple<shared_ptr<Vertex>, int>, bool>();

  for (auto ii = 0; ii < values.getSize(); ii++)
  {
    if (values[ii] >= 0.5)
    {
      int time = -1;
      for (auto &v : pool[ii]->path)
      {
        for (int delt = 0; delt <= instance->param.time_spacing; delt++)
          reserved[make_tuple(v, time + delt)] = true;
        time += 1;
      }
    }

    if (-PRIMAL_VALUE_EPS <= values[ii] && values[ii] <= PRIMAL_VALUE_EPS)
      continue;
    solution.push_back(make_tuple(pool[ii], values[ii]));
  }

  values.end();

  auto trt = make_shared<MasterResult>(obj_val,
                                       cardi_duals,
                                       solution);
  trt->num_col = cplex.getNcols();
  trt->num_row = cplex.getNrows();

  trt->i_at_cost_altered = i_at_cost_altered;
  trt->at_cost_altered = at_cost_altered;

  trt->reserved = reserved;

  if (instance->param.ctrs_string.find('V') != std::string::npos || instance->param.formul_string.find('V') != string::npos)
    trt->vertex_conflict_duals = vertex_conflict_duals;
  if (instance->param.ctrs_string.find('A') != std::string::npos || instance->param.formul_string.find('A') != string::npos)
    trt->arc_conflict_duals = arc_conflict_duals;
  if (instance->param.ctrs_string.find('T') != std::string::npos)
    trt->time_conflict_duals = time_conflict_duals;
  if (instance->param.ctrs_string.find('W') != std::string::npos || instance->param.formul_string.find('W') != string::npos)
    trt->wait_conflict_duals = wait_conflict_duals;
  if (instance->param.ctrs_string.find('R') != std::string::npos || instance->param.formul_string.find('R') != string::npos)
    trt->range_conflict_duals = range_conflict_duals;
  if (instance->param.ctrs_string.find('O') != std::string::npos)
    trt->old_conflict_duals = old_conflict_duals;
  if (instance->param.ctrs_string.find('C') != std::string::npos)
    trt->clique_conflict_duals = clique_conflict_duals;
  if (instance->param.ctrs_string.find('1') != std::string::npos)
    trt->tvc1_conflict_duals = tvc1_conflict_duals;
  if (instance->param.ctrs_string.find('2') != std::string::npos)
    trt->tvc2_conflict_duals = tvc2_conflict_duals;
  if (instance->param.ctrs_string.find('3') != std::string::npos)
    trt->tvc3_conflict_duals = tvc3_conflict_duals;
  if (instance->param.ctrs_string.find('4') != std::string::npos)
    trt->tvc4_conflict_duals = tvc4_conflict_duals;
  // Add conflict 9 sum up

  return trt;
}
