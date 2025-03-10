#ifndef SY_PROBLEM_DATA_H
#define SY_PROBLEM_DATA_H

#include "sy_unordered_map.h"
#include "sy_log.h"
#include "parameters.h"
#include "graph.h"

class Instance
{
public:
  // Instance & Algorithm parameters
  Param param;
  SYLog sylog;

  // Augmented directed graph
  Graph graph;
	
  Instance(Param param);

private:
  void read_scene_file(string file);
  void read_map_file(string mapfile);
  void log_instance_info();
};

#endif