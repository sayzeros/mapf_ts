#include "sub_solver.h"
#include "sy_log.h"
#include "spp_teg.h"

#include <set>
#include <string>

SubResult SubSolver::solve(int agent_id,
                           ColumnPool &node_pool,
                           shared_ptr<ColumnPool> global_pool,
                           shared_ptr<MasterResult> mas_result,
                           ErasedArcsOfAgents &erasedarcs)
{
  /* Solve SPP */
#ifdef RECORD_SPP_TIME
  auto spp_start = get_now();
#endif

  tuple<vector<shared_ptr<Vertex>>, double, double> result;
  if (instance->param.sub_algorithm == Param::SubAlg::dijkstra)
  {
    result = solve_spp_teg_dijkstra(instance->graph, agent_id, instance->param.num_agents, instance->param.time_spacing, mas_result, erasedarcs[agent_id]);
  }
  else if (instance->param.sub_algorithm == Param::SubAlg::astar)
  {
    result = solve_spp_teg_astar(instance->graph, agent_id, instance->param.num_agents, instance->param.time_spacing, mas_result, erasedarcs[agent_id]);
  }
  else
  {
    println("unimplemented");
    throw;
  }

#ifdef RECORD_SPP_TIME
  auto spp_time = get_elapsed(spp_start);
#ifdef SHOW_MASSUB_RESULT
  if (std::get<1>(result) < 0 - REDUCED_COST_EPS)
    println("            SP-{}: spp-time {:.3f}  best-reduced-cost {:.6f}", agent_id, spp_time, std::get<1>(result)); // 12 initial spaces
#endif
#endif
  /* Do somthing with the result */
  int num_new_col = 0;

  auto &[path, reduced_cost, path_cost] = result;
  if (reduced_cost < 0 - REDUCED_COST_EPS || global_pool->size() <= instance->param.num_agents) // profitable column founded
  {
    auto new_col = make_shared<Column>(instance, agent_id, path, path_cost);
    
    // set coefficents of separated cuts
    if (instance->param.ctrs_string.find('T') != std::string::npos)
    {
      for (auto [key, _]: mas_result->time_conflict_duals)
      {
        auto &[agent_id, v_id, time, delt] = key;
        double coeff = 0;
        if (new_col->agent_id == agent_id && time + 1 <= new_col->path.size()-1 && new_col->path[time+1]->id == v_id)
          coeff = 1;
        else if (new_col->agent_id != agent_id && time + delt + 1 <= new_col->path.size()-1 && new_col->path[time+delt+1]->id == v_id)
          coeff = 1;
        if (coeff > 0)
          new_col->time_conflict_coeff[key] = coeff;
      }
    }
    if (instance->param.ctrs_string.find('O') != std::string::npos)
    {
      for (auto [key, _]: mas_result->old_conflict_duals)
      {
        auto &[agent_id, v_id, time] = key;
        double coeff = 0;
        if (new_col->agent_id == agent_id && time + 1 <= new_col->path.size()-1 && new_col->path[time+1]->id == v_id)
          coeff = 1;
        else if (new_col->agent_id != agent_id)
        {
          for (int delt = 0; delt <= instance->param.time_spacing; delt++)
          {
            if (time + delt + 1 <= new_col->path.size()-1 && new_col->path[time+delt+1]->id == v_id)
              coeff += 1 / ((double) instance->param.time_spacing + 1);
          }
        }
        if (coeff > 0)
          new_col->old_conflict_coeff[key] = coeff;
      }
    }
    if (instance->param.ctrs_string.find('1') != std::string::npos)
    {
      for (auto [key, _]: mas_result->tvc1_conflict_duals)
      {
        auto &[agent_id, v_id, w_id, time] = key;
        double coeff = 0;
        if (new_col->agent_id == agent_id && time + 1 + 1 <= new_col->path.size()-1)
        {
          if (new_col->path[time+1]->id == v_id || new_col->path[time+1]->id == w_id || new_col->path[time+2]->id == v_id || new_col->path[time+2]->id == w_id)
          {
            coeff = 1;
          }
        }
        else if (new_col->agent_id != agent_id)
        {
          for (int delt = 0; delt <= instance->param.time_spacing - 1; delt++)
          {
            if (time + 1 + delt + 1 > new_col->path.size()-1)
              break;
            if ((new_col->path[time + 1 + delt]->id == v_id && new_col->path[time + 1 + delt + 1]->id == w_id) || (new_col->path[time + 1 + delt]->id == w_id && new_col->path[time + 1 + delt + 1]->id == v_id))
            {
              coeff = 1;
            }
          }
        }
        if (coeff > 0)
          new_col->tvc1_conflict_coeff[key] = coeff;
      }
    }
    if (instance->param.ctrs_string.find('2') != std::string::npos)
    {
      for (auto [key, _]: mas_result->tvc2_conflict_duals)
      {
        auto &[agent_id, v_id, w_id, time] = key;
        double coeff = 0;
        if (new_col->agent_id == agent_id && time + 1 + 1 <= new_col->path.size()-1)
        {
          if (new_col->path[time+1]->id == v_id || new_col->path[time+1]->id == w_id || new_col->path[time+2]->id == v_id || new_col->path[time+2]->id == w_id)
          {
            coeff = 1;
          }
        }
        else if (new_col->agent_id != agent_id)
        {
          for (int delt = 0; delt <= instance->param.time_spacing - 1; delt++)
          {
            if (time + 1 - delt < 0 || time + 1 - delt + 1 > new_col->path.size()-1)
              continue;
            if ((new_col->path[time + 1 - delt]->id == v_id && new_col->path[time + 1 - delt + 1]->id == w_id) || (new_col->path[time + 1 - delt]->id == w_id && new_col->path[time + 1 - delt + 1]->id == v_id))
            {
              coeff = 1;
            }
          }
        }
        if (coeff > 0)
          new_col->tvc2_conflict_coeff[key] = coeff;
      }
    }
    if (instance->param.ctrs_string.find('3') != std::string::npos)
    {
      for (auto [key, _]: mas_result->tvc3_conflict_duals)
      {
        auto &[agent_id, v_id, w_id, time] = key;
        double coeff = 0;
        if (new_col->agent_id == agent_id && time + 1 + 1 <= new_col->path.size()-1)
        {
          if (new_col->path[time + 1]->id == v_id || new_col->path[time + 2]->id == v_id || new_col->path[time + 2]->id == w_id)
            coeff += 1;
        }
        else if (new_col->agent_id != agent_id)
        {
          if (time + 2 <= new_col->path.size()-1 && new_col->path[time + 1]->id == v_id && new_col->path[time + 2]->id == w_id)
            coeff += 1;
          for (int delt = 1; delt <= instance->param.time_spacing - 1; delt++)
          {
            if (time + 1 + delt + 1 > new_col->path.size()-1)
              break;
            if ((new_col->path[time + 1 + delt]->id == v_id && new_col->path[time + 1 + delt + 1]->id == w_id) || (new_col->path[time + 1 + delt]->id == w_id && new_col->path[time + 1 + delt + 1]->id == v_id))
              coeff += 1;
          }
          if (time + 2 + instance->param.time_spacing <= new_col->path.size()-1 && new_col->path[time + 1 + instance->param.time_spacing]->id == v_id && new_col->path[time + 2 + instance->param.time_spacing]->id == w_id)
            coeff += 1;
        }
        if (coeff > 0)
          new_col->tvc3_conflict_coeff[key] = coeff;
      }
    }

    if (instance->param.ctrs_string.find('4') != std::string::npos)
    {
      for (auto [key, _]: mas_result->tvc4_conflict_duals)
      {
        auto &[agent_id, v_id, w_id, time] = key;
        double coeff = 0;
        if (new_col->agent_id == agent_id && time + 1 + 1 <= new_col->path.size()-1)
        {
          if (new_col->path[time + 1]->id == v_id || new_col->path[time + 1]->id == w_id || new_col->path[time + 2]->id == w_id)
            coeff += 1;
        }
        else if (new_col->agent_id != agent_id)
        {
          if (time + 2 <= new_col->path.size()-1 && new_col->path[time + 1]->id == v_id && new_col->path[time + 2]->id == w_id)
            coeff += 1;
          for (int delt = 1; delt <= instance->param.time_spacing - 1; delt++)
          {
            if (time + 1 - delt + 1 > new_col->path.size()-1 || time + 1 - delt < 0)
              continue;
            if ((new_col->path[time + 1 - delt]->id == v_id && new_col->path[time + 1 - delt + 1]->id == w_id) || (new_col->path[time + 1 - delt]->id == w_id && new_col->path[time + 1 - delt + 1]->id == v_id))
              coeff += 1;
          }
          if (time + 2 - instance->param.time_spacing <= new_col->path.size()-1 && time + 1 - instance->param.time_spacing >= 0 && new_col->path[time + 1 - instance->param.time_spacing]->id == v_id && new_col->path[time + 2 - instance->param.time_spacing]->id == w_id)
            coeff += 1;
        }
        if (coeff > 0)
          new_col->tvc4_conflict_coeff[key] = coeff;
      }
    }

#ifdef COLGEN_COLUMN_DUPLICATION_CHECK
    for (auto col : node_pool)
    {
      if (*col == *new_col)
      {
        println("****************************");
        println("****************************");
        println("*** Same Column DETECTED ***");
        println("****************************");
        println("****************************\n");
        println("reduced cost: {}  eps: {}\n", reduced_cost, REDUCED_COST_EPS);
        println("col 1 |{}", col->name());
        println("col 2 |{}", new_col->name());
        println("\n****************************");
        println("****************************");
        println("*** Same Column DETECTED ***");
        println("****************************");
        println("****************************");
        exit(1);
      }
    }
#endif
    node_pool.push_back(new_col);
    global_pool->push_back(new_col);
    num_new_col += 1;
  }

  return SubResult(num_new_col);
}
