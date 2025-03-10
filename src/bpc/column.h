#ifndef SY_COLUMN_H
#define SY_COLUMN_H

#include <memory>
#include <string>
#include <vector>

#include "instance.h"
#include "graph.h"

class Column
{
public:
  shared_ptr<Instance> instance;
  int agent_id;
  vector<shared_ptr<Vertex>> path; // single-agent plan

  double obj_coeff;                // cost of single-agent plan (path length)
  i_values cardi_coeff;
  vt_values vertex_conflict_coeff;
  at_values arc_conflict_coeff;
  ivtt_values time_conflict_coeff;
  vt_values wait_conflict_coeff;
  vt_values range_conflict_coeff;
  ivt_values old_conflict_coeff;
  sy_map<int, double> clique_conflict_coeff;
  ivwt_values tvc1_conflict_coeff;
  ivwt_values tvc2_conflict_coeff;
  ivwt_values tvc3_conflict_coeff;
  ivwt_values tvc4_conflict_coeff;
  // Add conflict 3 coef

  bool is_dummy;
  bool is_altered = false;

  Column() {}
  Column(shared_ptr<Instance> instance) : instance(instance) {}
  Column(shared_ptr<Instance> instance, int agent_id, vector<shared_ptr<Vertex>> path, double cost);

  string name();
  void make_dummy();

  bool operator==(const Column &other) const;
};

std::ostream &operator<<(std::ostream &out, const Column &c);

using ColumnPool = vector<shared_ptr<Column>>;

#endif