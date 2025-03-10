#ifndef mapfts_alg_tools_h_
#define mapfts_alg_tools_h_

#include "instance.h"

struct collision_t{
  int i;
  int j;
  int t_i;
  int t_j;
  shared_ptr<Vertex> v;
};

struct constraint_t{
  int i;
  int t1;
  int t2;
  shared_ptr<Vertex> v;
  bool operator<(const constraint_t& r) const
  {
    if (i != r.i) return i < r.i;
    if (t1 != r.t1) return t1 < r.t1;
    if (t2 != r.t2) return t2 < r.t2;
    return v < r.v;
  }
};

#endif