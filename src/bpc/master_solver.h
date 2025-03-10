#ifndef SY_MASTER_SOLVER_H
#define SY_MASTER_SOLVER_H

#include "instance.h"
#include "column.h"
#include "clique.h"

#include <memory>
#include <ilcplex/ilocplex.h>
#include "solution.h"
#include "master_result.h"

class MasterSolver
{
public:
  shared_ptr<Instance> instance;
  int num_added_columns;
  sy_map<tuple<int, int, int, int>, bool> time_conflict_is_needed;
  sy_map<tuple<int, int, int>, bool> old_conflict_is_needed;
  sy_map<tuple<int, int, int, int>, bool> tvc1_conflict_is_needed;
  sy_map<tuple<int, int, int, int>, bool> tvc2_conflict_is_needed;
  sy_map<tuple<int, int, int, int>, bool> tvc3_conflict_is_needed;
  sy_map<tuple<int, int, int, int>, bool> tvc4_conflict_is_needed;
  sy_map<tuple<int, int, int, int>, bool> time_conflict_is_added;

  MasterSolver() {}
  MasterSolver(shared_ptr<Instance> instance) : instance(instance) {}
  ~MasterSolver();

  shared_ptr<MasterResult> solve_lp(const ColumnPool &pool, string *model_dump = nullptr);
  void add_separated_cuts(const ColumnPool &pool);

  void reset_model();
  void generate_model(const ColumnPool &pool, bool solve_linear, string *model_dump);

  shared_ptr<vector<Clique>> cliques; // cliques for clique conflict constraints
private:
  bool is_just_reset = false;

  IloEnv env;
  IloModel model;
  IloNumVarArray vars;
  IloNumVarArray pat_vars;
  IloObjective obj;
  sy_map<int, IloRange> cardi_ctrs; // cardi[agent_id]

  sy_map<tuple<int, int>, IloRange> vertex_conflict_ctrs; // vc[(v, t)]
  sy_map<tuple<int, int>, bool> vertex_conflict_is_needed;
  sy_map<tuple<int, int>, bool> vertex_conflict_is_added;
  sy_map<tuple<int, int>, IloRange> arc_conflict_ctrs; // ac[(a,t)]
  sy_map<tuple<int, int>, bool> arc_conflict_is_needed;
  sy_map<tuple<int, int>, bool> arc_conflict_is_added;
  sy_map<tuple<int, int, int, int>, IloRange> time_conflict_ctrs; // tc[(i,v,t,delt)]
  sy_map<tuple<int, int>, IloRange> wait_conflict_ctrs; // wc[(v, t)]
  sy_map<tuple<int, int>, bool> wait_conflict_is_needed;
  sy_map<tuple<int, int>, bool> wait_conflict_is_added;
  sy_map<tuple<int, int>, IloRange> range_conflict_ctrs; // rc[(v, t)]
  sy_map<tuple<int, int>, bool> range_conflict_is_needed;
  sy_map<tuple<int, int>, bool> range_conflict_is_added;
  sy_map<tuple<int, int, int>, IloRange> old_conflict_ctrs; // oc[(v, t)]
  sy_map<tuple<int, int, int>, bool> old_conflict_is_added;
  sy_map<int, IloRange> clique_conflict_ctrs; // cc
  sy_map<int, bool> clique_conflict_is_needed;
  sy_map<int, bool> clique_conflict_is_added;
  sy_map<tuple<int, int, int, int>, IloRange> tvc1_conflict_ctrs; // 1c[(i,v,w,t)]
  sy_map<tuple<int, int, int, int>, bool> tvc1_conflict_is_added;
  sy_map<tuple<int, int, int, int>, IloRange> tvc2_conflict_ctrs; // 2c[(i,v,w,t)]
  sy_map<tuple<int, int, int, int>, bool> tvc2_conflict_is_added;
  sy_map<tuple<int, int, int, int>, IloRange> tvc3_conflict_ctrs; // 3c[(i,v,w,t)]
  sy_map<tuple<int, int, int, int>, bool> tvc3_conflict_is_added;
  sy_map<tuple<int, int, int, int>, IloRange> tvc4_conflict_ctrs; // 4c[(i,v,w,t)]
  sy_map<tuple<int, int, int, int>, bool> tvc4_conflict_is_added;
  // Add conflict 2 ctrs

  IloCplex cplex;

  using IloData = tuple<
      IloEnv,
      IloNumVarArray,
      sy_map<string, IloRange>,
      IloCplex>; // env, vars, balan_ctrs, cardi_ctrs, cplex
  void solve(const ColumnPool &pool, bool solve_linear, string *model_dump);
};

#endif