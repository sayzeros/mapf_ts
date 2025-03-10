#ifndef SY_CBS_T_H
#define SY_CBS_T_H

#include <vector>
#include <chrono>

#include "instance.h"
#include "sy_log.h"
#include "mapfts_alg_tools.h"


struct cbs_node_t{
  vector<vector<shared_ptr<Vertex>>> routes;
  cbs_node_t* parent;
  constraint_t ctr;
  double value;
  int idx;
  int depth;
  bool is_leaf = false;
};

struct cbs_node_t_cmp{
  bool operator()(cbs_node_t* t, cbs_node_t* u){
    return t->value > u->value;
  }
};

void cbs_ts(shared_ptr<Instance> instance);

#endif