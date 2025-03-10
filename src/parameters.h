#ifndef SY_PARAMETERS_H
#define SY_PARAMETERS_H

#include "sy_abbreviate.h"
#include "sy_util.h"

  ////////////////////////
 /// Macro Parameters ///
////////////////////////
#define MSG_LEVEL 2 // 0: final result;  1: B&B node info;  2: ColGen info;  3: Master/Subproblem info
// #define SHOW_BBNODE_RESULT_MINIMAL // B&B node info when bound updated

#define RECORD_SPP_TIME
#define DUMMY_COLUMN_COST  1e6
#define DUAL_VALUE_EPS     1e-9
#define PRIMAL_VALUE_EPS   1e-9
#define INTEGRALITY_EPS    1e-9
#define REDUCED_COST_EPS   1e-5
#define CUT_VIOLATION_EPS  1e-6
#define TRIVIAL_BOUND      1e9  // trivial bound that satisfies -this < lb and ub < this
#define NUM_THREAD 12 


#if (MSG_LEVEL == 0)
#define SHOW_MINIMAL_RESULT
#endif
#if (MSG_LEVEL >= 1)
#define SHOW_BBNODE_RESULT
#endif
#if (MSG_LEVEL >= 2)
#define SHOW_COLGEN_RESULT
#endif
#if (MSG_LEVEL >= 3)
#define SHOW_MASSUB_RESULT
#else
#undef RECORD_SPP_TIME
#endif

  ////////////////////////
 /// Type Definitions ///
////////////////////////
#define DOUBLE_NULL std::numeric_limits<double>::max()
typedef sy_map<int, double> i_values;  // idx: agnet_id
typedef sy_map<tuple<int, int>, double> vt_values;  // idx: vertex_id, timestep 
typedef sy_map<tuple<int, int, int>, double> ivt_values;  // idx: agent_id, vertex_id, timestep 
typedef sy_map<tuple<int, int, int, int>, double> ivwt_values;  // idx: agent_id, vertex_id, vertex_id, timestep 
typedef sy_map<tuple<int, int, int, int>, double> ivtt_values;  // idx: agent_id, vertex_id, timestep, delt 
typedef sy_map<tuple<int, int>, double> at_values;  // idx: arc_id, timestep 



/* Parameter class to store Instance & Algorithm parameters*/
class Param {
public:
  enum Alg {
    lp_cplex,
    mip_cplex,
    ccg,
    bpc,
    cbs,
    lbb,
    cnb,
  };
  static const string strings_alg[];

  enum SubAlg {
    dijkstra,
    astar
  };
  static const string strings_subalg[];

  enum PrimalHeur {
    global_cols, // solve MIP with global column pool
    selected_global_cols, // solve MIP with selected columns in global column pool (criteria: # included in node lp sol)
    local_cols, // solve MIP with local pool of the node 
    lp_sol_highest_cols, // choose highest lp solution value plan for each line
    lp_sol_cols // solve MIP with columns which have positive lp solution values
  };

  enum StepSizeRule {
    fixed1,
    fixed2,
    fixed3,
    fixed4,
    fixed5,
    dimini,
    polyak,
  };
  static const string strings_stepsize[];

  /* Instance file related */
  string instance_file;

  /* Instance parameters */
	int num_agents;
  int time_spacing;

  /* Model parameters */
  string ctrs_string;
  string formul_string;

  /* Algorithm parameters */
  double gap_limit = 0;	
	int time_limit = 60;
  int cut_rounds = 0; // if negative, no limit for the number of cut rounds
  int cut_depth = -1; // if negative, separate cuts in all nodes
	Alg algorithm;
  SubAlg sub_algorithm;
  StepSizeRule stepsize_rule;

  string alg_str;
  string subalg_str;
  string stepsize_str;
};

#endif