#ifndef SY_BNB_NODE_H
#define SY_BNB_NODE_H

#include "instance.h"
#include "column.h"
#include "master_solver.h"
#include "sub_solver.h"
// #include "cut_pool.h"
#include "branching_rule.h"

#include <memory>
#include <vector>

class BPCNode
{
public:
  shared_ptr<Instance> instance;
  shared_ptr<MasterResult> mas_result; // save master solver result to provide dual information for primal heuristic 

  shared_ptr<ColumnPool> pool;
  ColumnPool local_pool;

  /* Bound related */
  double lb;                                              // LP relaxation optimum obj val
  double ub;                                              // Upper bound from some feasible solution
  double parent_lb;                                       // LB of parent node, used to determine the priority in the queue.
  vector<tuple<shared_ptr<Column>, double>> base_columns; // optimal columns selected by the LP solver with the coefficients

  /* Branching related */
  ErasedArcsOfAgents erased_arcs;
  shared_ptr<BranchingRule> branching_rule;

  /* Node basic info */
  int depth;               // depth in the BnB tree
  int id;                  // node id
  static int num_explored; // # of BnBNodes explored
  int node_num_new_cols;
  int num_iter;

  /* #Cuts separated */
  int num_vc_ctrs;
  int num_ac_ctrs;
  int num_tc_ctrs;
  int num_wc_ctrs;
  int num_rc_ctrs;
  int num_oc_ctrs;
  int num_cc_ctrs;
  int num_1c_ctrs;
  int num_2c_ctrs;
  int num_3c_ctrs;
  int num_4c_ctrs;
  // Add conflict 0 #cut

  /* Time records */
  vector<double> all_times_spent_on_sub;
  double avg_time_spent_on_sp;
  double total_time_spent_on_sub;
  double total_time_spent_on_sepa;
  double total_time_spent_on_mas;
  double total_time_on_heur;
  double total_time_spent;
  double time_on_erasing_arcs;

  /* Constructor that will be called to generate the root node */
  BPCNode(shared_ptr<Instance> instance,
          shared_ptr<ColumnPool> pool,
          shared_ptr<MasterSolver> mas_solver, shared_ptr<SubSolver> sub_solver);

  /* Constructor to generate a child node from some parent node with some branching rule*/
  BPCNode(const BPCNode& parent,
          shared_ptr<BranchingRule> branching_rule,
          bool is_left_child);

  void solve();
  MAPFSol primal_heuristic(); 

  bool is_feasible();
  bool has_fractional_solution();

  void print_result(int num_unexplored, bool node_feasible, bool node_suboptimal, bool node_fractional, double global_ub);

private:
  shared_ptr<MasterSolver> mas_solver;
  shared_ptr<SubSolver> sub_solver;

  vector<tuple<int, int>> agent_order; // for primal heuristic

  void update_erased_arcs();
  void remove_incompatible_columns();

  shared_ptr<vector<Clique>> cliques; // cliques for clique conflict constraints
};

class BnBNodeCompare
{
public:
  bool operator()(const shared_ptr<BPCNode> &n1, const shared_ptr<BPCNode> &n2) const
  {
    return (n1->parent_lb > n2->parent_lb);
  }
};

#endif