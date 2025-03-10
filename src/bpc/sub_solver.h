#ifndef SY_SUB_SOLVER_H
#define SY_SUB_SOLVER_H

#include "instance.h"
#include "column.h"
#include "master_solver.h"

#include <memory>
#include <map>
#include <vector>

struct SubResult
{
  int num_new_cols;
  bool is_infeasible = false;
  SubResult(int num_new_cols) : num_new_cols(num_new_cols){};
};

class SubSolver
{
public:
  shared_ptr<Instance> instance;

  SubSolver() {}
  SubSolver(shared_ptr<Instance> instance) : instance(instance){};

  SubResult solve(int agent_id,
                  ColumnPool &node_pool,
                  shared_ptr<ColumnPool> global_pool,
                  shared_ptr<MasterResult> mas_result,
                  ErasedArcsOfAgents &erasedarcs);

private:
};

#endif