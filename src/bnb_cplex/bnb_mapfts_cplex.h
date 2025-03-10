#ifndef _BNB_CPLEX_H_
#define _BNB_CPLEX_H_

#include "instance.h"
#include "sy_log.h"


class BnBSolver {
public: 
  shared_ptr<Instance> instance;
  
  BnBSolver() {}
  BnBSolver(shared_ptr<Instance> instance) : instance(instance) {}

  void solve(bool solve_linear=true, bool write_solution=false);
  void write_solution();
};

#endif