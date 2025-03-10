#include "solution.h"
#include "sy_log.h"
#include <fmt/os.h>


void MAPFSol::print()
{
  fmt::print("obj: {}\n", obj_val);
  int i = 0;
  for (auto &path : paths)
  {
    fmt::print("Agent {}, cost {}, path ", i, (int) path.size() - 1);
    int cnt = 0;
    for (auto &coord : path)
    {
      fmt::print("{}", coord.name());
      if (cnt >= path.size() - 1)
      {
        break;
      }
      fmt::print(",");
      cnt += 1;
    }
    fmt::print("\n");
    i += 1;
  }
}

string MAPFSol::write()
{
  string filename = fmt::format(
      "{}-{}agents-{}ts-{}ctrs-{}.sol",
      inst->param.instance_file.substr(0, inst->param.instance_file.find_last_of('.')),
      inst->param.num_agents,
      inst->param.time_spacing,
      inst->param.ctrs_string,
      inst->param.alg_str);
  string output_path;
  if (filename[0] == '/')
  {
    string soldir = filename.substr(0, filename.find_last_of('/'));
    soldir = soldir.substr(0, soldir.find_last_of('/'));
    soldir = soldir.substr(0, soldir.find_last_of('/'));
    string filepath = filename.substr(0, filename.find_last_of('/'));
    filename = filepath.substr(filepath.find_last_of('/') + 1) + "/" + filename.substr(filename.find_last_of('/') + 1);
    output_path = soldir + "/output/solution/" + filename;
  }
  else
  {
    output_path = "./output/solution/" + filename;
  }
  auto out = fmt::output_file(output_path);
  out.print("{}\n\n", obj_val);
  int i = 0;
  for (auto &path : paths)
  {
    out.print("Agent {}, cost {}, path ", i, path.size() - 1);
    int cnt = 0;
    for (auto &coord : path)
    {
      out.print("{}", coord.name());
      if (cnt >= path.size() - 1)
      {
        break;
      }
      out.print(",");
      cnt += 1;
    }
    out.print("\n");
    i += 1;
  }
  return output_path;
}

VertexConflict MAPFSol::find_first_vertex_conflict()
{
  int num_agents = paths.size();
  int makespan = 0;
  for (auto &path : paths)
    makespan = makespan < path.size() - 1 ? path.size() - 1 : makespan;

  sy_map<tuple<int, int, int>, int> vt_visited_by; // input: (x,y,t) / output: i
  for (int t = 0; t <= makespan; t++)
  {
    for (int i = 0; i < num_agents; i++)
    {
      if (t < paths[i].size())
      {
        if (vt_visited_by[make_tuple(paths[i][t].x, paths[i][t].y, t)] == 0)
          vt_visited_by[make_tuple(paths[i][t].x, paths[i][t].y, t)] = i + 1;
        else
          return VertexConflict(vt_visited_by[make_tuple(paths[i][t].x, paths[i][t].y, t)] - 1, i, paths[i][t].x, paths[i][t].y, t);
      }
    }
  }
  return VertexConflict();
}

ArcConflict MAPFSol::find_first_arc_conflict()
{
  int num_agents = paths.size();
  int makespan = 0;
  for (auto &path : paths)
    makespan = makespan < path.size() - 1 ? path.size() - 1 : makespan;

  sy_map<tuple<int, int, int, int, int>, int> vwt_visited_by; // input: (v_x,v_y,w_x,w_y,t) / output: i
  for (int t = 0; t <= makespan; t++)
  {
    for (int i = 0; i < num_agents; i++)
    {
      if (t < paths[i].size() - 1)
      {
        if (vwt_visited_by[make_tuple(paths[i][t].x, paths[i][t].y, paths[i][t + 1].x, paths[i][t + 1].y, t)] == 0)
          vwt_visited_by[make_tuple(paths[i][t].x, paths[i][t].y, paths[i][t + 1].x, paths[i][t + 1].y, t)] = i + 1;
        else
          return ArcConflict(vwt_visited_by[make_tuple(paths[i][t].x, paths[i][t].y, paths[i][t + 1].x, paths[i][t + 1].y, t)] - 1, i, paths[i][t].x, paths[i][t].y, paths[i][t + 1].x, paths[i][t + 1].y, t);
        if (paths[i][t].x != paths[i][t + 1].x || paths[i][t].y != paths[i][t + 1].y)
        {
          if (vwt_visited_by[make_tuple(paths[i][t + 1].x, paths[i][t + 1].y, paths[i][t].x, paths[i][t].y, t)] == 0)
            vwt_visited_by[make_tuple(paths[i][t + 1].x, paths[i][t + 1].y, paths[i][t].x, paths[i][t].y, t)] = i + 1;
          else
            return ArcConflict(vwt_visited_by[make_tuple(paths[i][t + 1].x, paths[i][t + 1].y, paths[i][t].x, paths[i][t].y, t)] - 1, i, paths[i][t + 1].x, paths[i][t + 1].y, paths[i][t].x, paths[i][t].y, t);
        }
      }
    }
  }
  return ArcConflict();
}

TimeConflict MAPFSol::find_first_time_conflict()
{
  int &ts = inst->param.time_spacing;
  int num_agents = paths.size();
  int makespan = 0;
  for (auto &path : paths)
    makespan = makespan < path.size() - 1 ? path.size() - 1 : makespan;

  sy_map<tuple<int, int, int>, tuple<int, int>> vt_occupied_by; // input: (x,y,t) / output: (i, t_i)
  for (int t = 0; t <= makespan; t++)
  {
    for (int i = 0; i < num_agents; i++)
    {
      if (t < paths[i].size())
        for (int tt = t; tt <= makespan && tt <= t + ts; tt++)
        {
          if (std::get<0>(vt_occupied_by[make_tuple(paths[i][t].x, paths[i][t].y, tt)]) == 0)
            vt_occupied_by[make_tuple(paths[i][t].x, paths[i][t].y, tt)] = make_tuple(i + 1, t);
          else if (std::get<0>(vt_occupied_by[make_tuple(paths[i][t].x, paths[i][t].y, tt)]) == i + 1)
            vt_occupied_by[make_tuple(paths[i][t].x, paths[i][t].y, tt)] = make_tuple(i + 1, t);
          else
            return TimeConflict(std::get<0>(vt_occupied_by[make_tuple(paths[i][t].x, paths[i][t].y, tt)]) - 1, i, paths[i][t].x, paths[i][t].y, std::get<1>(vt_occupied_by[make_tuple(paths[i][t].x, paths[i][t].y, tt)]), t);
        }
    }
  }
  return TimeConflict();
}