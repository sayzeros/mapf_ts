#ifndef SY_MASTER_RESULT_H
#define SY_MASTER_RESULT_H

#include "instance.h"

struct MasterResult
{
  double obj_val;
  i_values cardi_duals;

  vt_values vertex_conflict_duals;
  at_values arc_conflict_duals;
  ivtt_values time_conflict_duals;
  vt_values wait_conflict_duals;
  vt_values range_conflict_duals;
  ivt_values old_conflict_duals;
  sy_map<int, double> clique_conflict_duals;
  ivwt_values tvc1_conflict_duals;
  ivwt_values tvc2_conflict_duals;
  ivwt_values tvc3_conflict_duals;
  ivwt_values tvc4_conflict_duals;
  // Add conflict 1 decl

  sy_map<tuple<shared_ptr<Vertex>, int>, bool> reserved;

  vector<tuple<shared_ptr<Column>, double>> solution;

  int num_row;
  int num_col;

  sy_map<int, sy_map<tuple<shared_ptr<Arc>, int>, double>> i_at_cost_altered;
  sy_map<tuple<shared_ptr<Arc>, int>, double> at_cost_altered;

  MasterResult() {}

  // no ctrs
  MasterResult(double obj_val,
               i_values &cardi_duals,
               vector<tuple<shared_ptr<Column>, double>> solution)
      : obj_val(obj_val),
        cardi_duals(cardi_duals),
        solution(solution)
  {
  }

  void print_solution()
  {
    println("{}", string(60, '='));
    println("LP Solution");
    println("{}", string(60, '-'));
    println("obj : {}", obj_val);
    for (auto [col, val] : solution)
    {
      if (val > 0)
      {
        println("val {} {}", val, col->name());
      }
    }
    println("{}", string(60, '='));
  }
};

#endif