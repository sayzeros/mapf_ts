#ifndef SY_CUT_POOL_H
#define SY_CUT_POOL_H

#include "instance.cpp"

struct CutPool {
  vector<tuple<int, int, int, int>> tc_keys; // agent_id, vertex_id, time, delta_t
};

#endif 
