#include "bnb_mapfts_cplex.h"
#include "solution.h"

#include <ilcplex/ilocplex.h>

void BnBSolver::solve(bool solve_linear, bool write_solution)
{
  time_point mip_start_time = get_now();
  std::stringstream temp_name; // to give name for columns and rows

  auto &graph = instance->graph;
  auto &sylog = instance->sylog;
  auto &ts = instance->param.time_spacing;

  int num_agents = instance->param.num_agents;
  int len_horizon = 16;
  int width = graph.width;
  int height = graph.height;
  int num_arcs = graph.arcs.size();

  IloEnv env;
  IloModel model(env);
  IloNumVarArray vars(env);
  IloExpr expr;

  /* Generate decision variables (row-wise) */
  debugln("Generating decision variables");
  sy_map<int, sy_map<int, sy_map<int, IloNumVar>>> xx; 

  for (int i = 0; i < num_agents; i++)
  {
    for (int t = 0; t < len_horizon; t++)
    {
      for (int a_id = 0; a_id < num_arcs; a_id++)
      {
        if (graph.get_arc(a_id)->from->type == Vertex::Type::DUMMY_ORIGIN)
        {
          if (t != 0 || i != graph.get_arc(a_id)->from->agent_id)
          {
            temp_name.str("");
            temp_name << fmt::format("x_{}_{}_{}", i, graph.get_arc(a_id)->name(), t);
            xx[i][a_id][t] =
                IloNumVar(env, 0, 0,
                          solve_linear ? IloNumVar::Float : IloNumVar::Int, temp_name.str().c_str());
            vars.add(xx[i][a_id][t]);
            continue;
          }
        }
        if (graph.get_arc(a_id)->to->type == Vertex::Type::DUMMY_DESTINATION)
        {
          if (i != graph.get_arc(a_id)->to->agent_id)
          {
            temp_name.str("");
            temp_name << fmt::format("x_{}_{}_{}", i, graph.get_arc(a_id)->name(), t);
            xx[i][a_id][t] =
                IloNumVar(env, 0, 0,
                          solve_linear ? IloNumVar::Float : IloNumVar::Int, temp_name.str().c_str());
            vars.add(xx[i][a_id][t]);
            continue;
          }
        }
        temp_name.str("");
        temp_name << fmt::format("x_{}_{}_{}", i, graph.get_arc(a_id)->name(), t);
        xx[i][a_id][t] =
            IloNumVar(env, 0, 1,
                      solve_linear ? IloNumVar::Float : IloNumVar::Int, temp_name.str().c_str());
        vars.add(xx[i][a_id][t]);
      }
    }
  }
  /* Declare objective */
  debugln("Generating objective function");
  expr = IloExpr(env, 0);
  for (int i = 0; i < num_agents; i++)
  {
    for (int t = 0; t < len_horizon; t++)
    {
      for (int a_id = 0; a_id < num_arcs; a_id++)
      {
        expr += xx[i][a_id][t] * graph.get_arc(a_id)->cost;
      }
    }
  }
  IloObjective obj = IloMinimize(env, expr);

  /* Declare constraints */
  sy_map<int, sy_map<int, sy_map<int, IloRange>>> fc_ctrs; // fc[v][t][i]
  sy_map<int, IloRange> fc_o_ctrs;                         // fc_o[i]
  sy_map<int, IloRange> fc_d_ctrs;                         // fc_d[i]

  /* Generate constraints (only with rhs constants; lhs columns will be added later) */
  IloRangeArray ctrs(env);
  int num_ctrs = 0;

  debugln("Generating fc constraints");
  // fc constraints (for original vertices)
  for (int x = 0; x < width; x++)
  {
    for (int y = 0; y < height; y++)
    {
      auto v = graph.get_vertex(Coord(x, y));
      if (v == nullptr)
        continue;
      for (int t = 0; t < len_horizon; t++)
      {
        for (int i = 0; i < num_agents; i++)
        {
          expr = IloExpr(env, 0);
          temp_name.str("");
          temp_name << fmt::format("fc_{}_{}_{}", v->name(), t, i);
          for (auto &a : graph.get_out_arcs(v))
          {
            expr += xx[i][a->id][t];
          }
          for (auto &a : graph.get_in_arcs(v))
          {
            if (t - a->cost >= 0)
              expr -= xx[i][a->id][t - a->cost];
          }
          fc_ctrs[v->id][t][i] = IloRange(env, 0, expr, 0, temp_name.str().c_str());
          ctrs.add(fc_ctrs[v->id][t][i]);
          num_ctrs += 1;
        }
      }
    }
  }

  for (int i = 0; i < num_agents; i++)
  {
    expr = IloExpr(env, 0);
    temp_name.str("");
    temp_name << fmt::format("fco_{}", i);
    for (auto &a : graph.get_out_arcs(graph.get_dummy_vertex(Vertex::Type::DUMMY_ORIGIN, i)))
    {
      expr += xx[i][a->id][0];
    }
    fc_o_ctrs[i] = IloRange(env, 1, expr, 1, temp_name.str().c_str());
    ctrs.add(fc_o_ctrs[i]);
    num_ctrs += 1;
  }

  for (int i = 0; i < num_agents; i++)
  {
    expr = IloExpr(env, 0);
    temp_name.str("");
    temp_name << fmt::format("fcd_{}", i);
    for (auto &a : graph.get_in_arcs(graph.get_dummy_vertex(Vertex::Type::DUMMY_DESTINATION, i)))
    {
      for (int t = 0; t < len_horizon; t++)
      {
        expr += xx[i][a->id][t - a->cost]; // Note a->cost == 0
      }
    }
    fc_d_ctrs[i] = IloRange(env, 1, expr, 1, temp_name.str().c_str());
    ctrs.add(fc_d_ctrs[i]);
    num_ctrs += 1;
  }

  /* Conflict constraints (Constraints) */ 
  // Vertex constraints
  if (instance->param.ctrs_string.find('V') != string::npos)
  {
    for (int x = 0; x < width; x++)
    {
      for (int y = 0; y < height; y++)
      {
        auto v = graph.get_vertex(Coord(x, y));
        if (v == nullptr)
          continue;
        for (int t = 0; t < len_horizon; t++)
        {
          expr = IloExpr(env, 0);
          temp_name.str("");
          temp_name << fmt::format("vc_{}_{}", v->name(), t);
          for (int i = 0; i < num_agents; i++)
          {
            for (auto a : graph.get_out_arcs(v))
              expr += xx[i][a->id][t];
          }
          ctrs.add(IloRange(env, -1, -expr));
          ctrs[num_ctrs].setName(temp_name.str().c_str());
          num_ctrs += 1;
        }
      }
    }
  }
  // Arc constraints
  if (instance->param.ctrs_string.find('A') != string::npos)
  {
    for (int x = 0; x < width; x++)
    {
      for (int y = 0; y < height; y++)
      {
        auto v = graph.get_vertex(Coord(x, y));
        if (v == nullptr)
          continue;
        for (auto a : graph.get_outout_arcs(v))
        {
          auto w = a->to;
          if (w->type == Vertex::Type::DUMMY_DESTINATION)
            continue;
          for (int t = 0; t < len_horizon; t++)
          {
            expr = IloExpr(env, 0);
            temp_name.str("");
            temp_name << fmt::format("ac_{}_{}_{}", v->name(), w->name(), t);
            for (int i = 0; i < num_agents; i++)
            {
              expr += xx[i][a->id][t];
              expr += xx[i][graph.get_arc(w, v)->id][t];
            }
            ctrs.add(IloRange(env, -1, -expr));
            ctrs[num_ctrs].setName(temp_name.str().c_str());
            num_ctrs += 1;
          }
        }
      }
    }
  }
  // Old-time-spacing constraints
  if (instance->param.ctrs_string.find('O') != string::npos)
  {
    for (int x = 0; x < width; x++)
    {
      for (int y = 0; y < height; y++)
      {
        auto v = graph.get_vertex(Coord(x, y));
        if (v == nullptr)
          continue;
        for (int t = 0; t < len_horizon; t++)
        {
          for (int i = 0; i < num_agents; i++)
          {
            expr = IloExpr(env, 0);
            temp_name.str("");
            temp_name << fmt::format("otc_{}_{}_{}", v->name(), i, t);
            for (auto a : graph.get_out_arcs(v))
            {
              expr += xx[i][a->id][t] * (ts + 1);
              for (int j = 0; j < num_agents; j++)
              {
                if (j == i)
                  continue;
                for (int tt = t; tt < len_horizon && tt <= t + ts; tt++)
                  expr += xx[j][a->id][tt];
              }
            }
            ctrs.add(IloRange(env, -ts - 1, -expr));
            ctrs[num_ctrs].setName(temp_name.str().c_str());
            num_ctrs += 1;
          }
        }
      }
    }
  }
  // Time-spacing constraints
  if (instance->param.ctrs_string.find('T') != string::npos)
  {
    for (int x = 0; x < width; x++)
    {
      for (int y = 0; y < height; y++)
      {
        auto v = graph.get_vertex(Coord(x, y));
        if (v == nullptr)
          continue;
        for (int t = 0; t < len_horizon; t++)
        {
          for (int delt = 0; delt + t < len_horizon && delt <= ts; delt++)
            for (int i = 0; i < num_agents; i++)
            {
              expr = IloExpr(env, 0);
              temp_name.str("");
              temp_name << fmt::format("tc_{}_{}_{}_{}", v->name(), i, t, delt);
              for (auto a : graph.get_out_arcs(v))
              {
                expr += xx[i][a->id][t];
                for (int j = 0; j < num_agents; j++)
                {
                  if (j == i)
                    continue;
                  expr += xx[j][a->id][t + delt];
                }
              }
              ctrs.add(IloRange(env, -1, -expr));
              ctrs[num_ctrs].setName(temp_name.str().c_str());
              num_ctrs += 1;
            }
        }
      }
    }
  }

  /* Generate model and solve */
  debugln("Add all to the model; {}", get_elapsed(mip_start_time));
  // objective
  model.add(obj);
  // constraints
  model.add(ctrs);

  /* Make cplex! */
  IloCplex cplex(model);


  /* Save .lp file */
  string output_name = "test.lp";
  try
  {
    cplex.exportModel(output_name.c_str());
  }
  catch (IloException &e)
  {
    std::cerr << "Export IloException: " << e << std::endl;
  }

  cplex.setParam(IloCplex::Param::TimeLimit, instance->param.time_limit);
#ifdef DEBUG
  cplex.setParam(IloCplex::Param::MIP::Display, 3);
#else
  cplex.setOut(env.getNullStream());
#endif

  time_point mip_load_time = get_now();
  sylog.add_alg_result("time_load", get_elapsed(mip_start_time, mip_load_time));

  debugln("Start to solve!!");
  debugnewline(2);

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

  debugnewline(2);
  debugln("Done solving; {}", get_elapsed(mip_start_time));

  /* Interpret solver results */
  double obj_val;
  if (!solved)
  {
    debugln(">> Integer feasible solution not found (cplex)");
    obj_val = DOUBLE_NULL;
  }
  else
  {
    obj_val = cplex.getObjValue();
  }
  string prob_status;
  switch (cplex.getStatus())
  {
  case 0:
    prob_status = "unknown";
    break;
  case 1:
    prob_status = "feasible";
    break;
  case 2:
    prob_status = "optimal";
    break;
  case 3:
    prob_status = "infeasible";
    break;
  case 4:
    prob_status = "unbounded";
    break;
  case 5:
    prob_status = "inf_or_ubd";
    break;
  case 6:
    prob_status = "eror";
    break;
  default:
    prob_status = "idk";
    break;
  }
  if (solve_linear)
  {
    sylog.add_alg_result("obj", "-");
    sylog.add_alg_result("ub", "-");
    sylog.add_alg_result("lb", obj_val);
    sylog.add_alg_result("#node", "-");
  }
  else
  {
    sylog.add_alg_result("obj", obj_val);
    sylog.add_alg_result("ub", obj_val);
    sylog.add_alg_result("lb", cplex.getBestObjValue());
    sylog.add_alg_result("#node", cplex.getNnodes());
  }
  sylog.add_alg_result("time_solve", get_elapsed(mip_load_time));
  sylog.add_alg_result("time", get_elapsed(mip_start_time));
  sylog.add_alg_result("stat", prob_status);
  sylog.add_alg_result("#row", cplex.getNcols());
  sylog.add_alg_result("#col", cplex.getNrows());
  debug_assert(vars.getSize() == cplex.getNcols());
  debug_assert(num_ctrs == cplex.getNrows());

  sylog.print_result();


  if (instance->param.alg_str == "mip_cplex")
  {
    /* Convert solution as MasRes format */
    if (write_solution && (prob_status == "optimal" || prob_status == "feasible"))
    { 
      double value;
      vector<vector<Coord>> paths;
      for (int i = 0; i < num_agents; i++)
      {
        double path_cost = 0;
        double dummy_origin_cnt = 0;
        double dummy_destination_cnt = 0;
        vector<Coord> v_path;
        for (int t = 0; t < len_horizon; t++)
        {
          for (int a_id = 0; a_id < num_arcs; a_id++)
          {
            value = cplex.getValue(xx[i][a_id][t]);
            if (value > 0)
            {
              if (graph.get_arc(a_id)->from->type == Vertex::Type::DUMMY_ORIGIN)
              {
                dummy_origin_cnt += value;
              }
              else
              {
                v_path.push_back(graph.get_arc(a_id)->from->coord);
              }
              if (graph.get_arc(a_id)->to->type == Vertex::Type::DUMMY_DESTINATION)
              {
                dummy_destination_cnt += value;
              }
              debugln("x_{}_{}_{}: {}", i, graph.get_arc(a_id)->name(), t, value);
              path_cost += value * graph.get_arc(a_id)->cost;
            }
          }
        }
        debug_assert(dummy_origin_cnt - 1 <= 1e-6);
        debug_assert(dummy_destination_cnt - 1 <= 1e-6);
        paths.push_back(v_path);
      }
      MAPFSol sol(instance, obj_val, paths);
      string write_path = sol.write();
      debugln("Save solution at {}", write_path);
      // check conflicts
      auto vc = sol.find_first_vertex_conflict();
      println("vc {}", vc.name());

      auto ac = sol.find_first_arc_conflict();
      println("ac {}", ac.name());

      auto tc = sol.find_first_time_conflict();
      println("tc {}", tc.name());
    }
  }

  env.end();
}
