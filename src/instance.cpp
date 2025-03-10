#include "instance.h"
#include "sy_util.h"

#include <fstream>

string get_str_token(string &line, string &delim);

Instance::Instance(Param param) : param(param)
{
  read_scene_file(param.instance_file);
  log_instance_info();
}

void Instance::log_instance_info()
{
  string &scenefile = param.instance_file;
  sylog.add_inst_info("inst_file", scenefile);
  string map = scenefile.substr(scenefile.find_last_of("/") + 1);
  map = map.substr(0, map.find_last_of("-"));
  map = map.substr(0, map.find_last_of("-"));
  sylog.add_inst_info("map", map);
  string type_id = scenefile.substr(scenefile.substr(0, scenefile.find_last_of("-")).find_last_of("-") + 1);
  type_id = type_id.substr(0, type_id.find_last_of("."));
  sylog.add_inst_info("type", type_id.substr(0, type_id.find_last_of("-")));
  sylog.add_inst_info("id", type_id.substr(type_id.find_last_of("-") + 1));

  sylog.add_inst_info("#agents", param.num_agents);
  sylog.add_inst_info("t_s", param.time_spacing);
  sylog.add_inst_info("formul", param.formul_string);
  sylog.add_inst_info("ctrs", param.ctrs_string);
  sylog.add_inst_info("time_lim", param.time_limit);
  sylog.add_inst_info("gap_lim", param.gap_limit);
  sylog.add_inst_info("alg", param.alg_str);
  sylog.add_inst_info("sub_alg", param.subalg_str);
  sylog.add_inst_info("stepsize_rule", param.stepsize_str);

  sylog.add_inst_info("#vert", to_string(graph.vertices.size()));
  sylog.add_inst_info("#vert_real", to_string(graph.vertices.size() - 2 * param.num_agents));
  sylog.add_inst_info("#vert_dumm", to_string(2 * param.num_agents));
  sylog.add_inst_info("#arc", to_string(graph.arcs.size()));
  sylog.add_inst_info("#arc_real", to_string(graph.arcs.size() - 2 * param.num_agents));
  sylog.add_inst_info("#arc_dumm", to_string(2 * param.num_agents));
}

void Instance::read_scene_file(string scenefile)
{
  string instance_foler = "./data/";
  std::ifstream instance_scene;
  if (scenefile[0] == '/')
  {
    instance_scene = std::ifstream(scenefile);
  }
  else
  {
    instance_scene = std::ifstream(instance_foler + scenefile);
  }

  string line;
  if (instance_scene.is_open())
  {
    getline(instance_scene, line);
    // D(cout << "\t" << line << endl); // First line should be version info
    int agent_id = 0;
    bool map_read = false;
    while (getline(instance_scene, line)) // iteration till eof
    {
      // D(cout << "\t" << line << endl);
      string delim = "\t";
      string token;
      Coord o, d;
      token = get_str_token(line, delim); // token1:  ?
      token = get_str_token(line, delim); // token2: mapfile name
      if (!map_read)
      {
        string mapfile = scenefile.substr(0, scenefile.find_last_of("/")) + "/" + token;
        if (mapfile[0] == '/')
        {
          read_map_file(mapfile);
        }
        else
        {
          read_map_file(instance_foler + mapfile);
        }
        map_read = true;
      }
      token = get_str_token(line, delim); // token3: width of the map
      token = get_str_token(line, delim); // token4: heights of the map
      token = get_str_token(line, delim); // token5: origin x-coordination
      o.x = stoi(token);
      token = get_str_token(line, delim); // token6: origin y-coordination
      o.y = stoi(token);
      token = get_str_token(line, delim); // token7: destination x-coordination
      d.x = stoi(token);
      token = get_str_token(line, delim); // token8: destination y-coordination
      d.y = stoi(token);
      token = get_str_token(line, delim); // token9: ?

      graph.add_dummy_vertex(Vertex::Type::DUMMY_ORIGIN, o, agent_id);
      graph.add_dummy_vertex(Vertex::Type::DUMMY_DESTINATION, d, agent_id);
      agent_id += 1;

      if (agent_id >= param.num_agents)
        break;
    }
  }
}

void Instance::read_map_file(string mapfile)
{
  std::ifstream map_file(mapfile);
  string line, token, delim;
  delim = " ";

  if (map_file.is_open())
  {
    getline(map_file, line); // First line should be type of the map. We ignore that information.
    getline(map_file, line); // First token of second line is "height"
    token = get_str_token(line, delim);
    graph.height = stoi(get_str_token(line, delim));
    getline(map_file, line); // First token of third line is "width"
    token = get_str_token(line, delim);
    graph.width = stoi(get_str_token(line, delim));
    getline(map_file, line); // Forth line is "map"

    for (int y = 0; y < graph.height; y++)
    {
      getline(map_file, line);
      for (int x = 0; x < graph.width; x++)
      {
        // Note: "." or "S" denote passable vertex. Otherwise, impassable
        if (line.at(x) == '.' or line.at(x) == 'S')
        {
          auto v = graph.add_vertex(Coord(x, y));
          graph.add_arc(v, v);
          auto left = graph.get_vertex(Coord(x - 1, y));
          auto down = graph.get_vertex(Coord(x, y - 1));
          if (x - 1 >= 0 && left != nullptr)
          {
            graph.add_arc(left, v);
            graph.add_arc(v, left);
          }
          if (y - 1 >= 0 && down != nullptr)
          {
            graph.add_arc(down, v);
            graph.add_arc(v, down);
          }
        }
        else
        {
          continue;
        }
      }
    }
  }
  else
  {
    println("Couldn't open the map file: maybe not exist");
  }
  map_file.close();
}

string get_str_token(string &line, string &delim)
{
  string token;
  int del_pos = line.find(delim);
  if (del_pos == -1)
  {
    token = line;
    line = "";
  }
  else
  {
    token = line.substr(0, del_pos);
    line = line.substr(del_pos + delim.length(), -1);
  }
  return token;
}