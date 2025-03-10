#ifndef SY_SOLUTION_H
#define SY_SOLUTION_H

#include "instance.h"
#include "conflict.h"
#include "column.h"

class MAPFSol
{
public:
  shared_ptr<Instance> inst;
  double obj_val;
  vector<vector<Coord>> paths;
  int ts;

  MAPFSol(shared_ptr<Instance> inst,
          double obj_val,
          vector<vector<Coord>> paths)
      : inst(inst), obj_val(obj_val), paths(paths) {}

  void print();
  string write();

  VertexConflict find_first_vertex_conflict();
  ArcConflict find_first_arc_conflict();
  TimeConflict find_first_time_conflict();
};

#endif