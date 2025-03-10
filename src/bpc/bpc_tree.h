#ifndef SY_BNB_TREE_H
#define SY_BNB_TREE_H

#include <memory>
#include <queue>
#include <vector>
#include <chrono>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "instance.h"
#include "bpc_node.h"
#include "sy_log.h"

class BPCTree
{
  using NodeQueue = std::priority_queue<
      shared_ptr<BPCNode>,
      vector<shared_ptr<BPCNode>>,
      BnBNodeCompare>;

  enum Status
  {
    ONGOING,
    TIME_LIMIT_REACHED,
    OPTIMAL_FOUNDED,
    INFEASIBLE,
    LP_OPTIMAL,
  };

public:
  shared_ptr<Instance> instance;
  shared_ptr<ColumnPool> pool;
  NodeQueue unexplored_nodes;
  Status status;

  /* Statistics maintained by update_stats */
  double ub;
  shared_ptr<BPCNode> node_attaining_ub;
  double root_lb;
  double gap_at_root;
  double gap;
  int max_depth;
  double elapsed_time;

  double total_time_on_mp;
  double total_time_on_sp;
  double total_time_on_sepa;
  double total_time_on_branching;
  double total_time_on_heur;
  double total_time_on_node;

  double avg_time_on_sp;
  double avg_time_on_erasing;
  int num_spp;

#ifdef SHOW_BBNODE_RESULT_MINIMAL
  double last_lb;
  double last_ub;
#endif

  // #Cuts separated
  int num_vc_ctrs, num_ac_ctrs;
  int num_tc_ctrs;
  int num_wc_ctrs;
  int num_rc_ctrs;
  int num_oc_ctrs;
  int num_cc_ctrs;
  int num_1c_ctrs;
  int num_2c_ctrs;
  int num_3c_ctrs;
  int num_4c_ctrs;

  /* Statistics updated manually */
  int num_nodes_generated; // update when generate node (construct root & branching)
  double lb;               // (after solving root & retrieve from the queue)

  /* Constructor */
  BPCTree() {}
  BPCTree(shared_ptr<Instance> instance);

  /* Public methods */
  void solve_root(bool called_by_user = true);
  void explore_tree();

private:
  shared_ptr<MasterSolver> mas_solver;
  shared_ptr<SubSolver> sub_solver;

  bool bound_updated;
  time_point tree_start_time;

  void update_stats(shared_ptr<BPCNode> node);

  void branch(shared_ptr<BPCNode> node);
  bool branch_on_forbid_arcs(shared_ptr<BPCNode> node);

  void primal_heuristic(shared_ptr<BPCNode> current_node);

  void print_result();

  string get_status_string();
};

#endif