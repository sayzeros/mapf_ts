#include "sy_log.h"

#include "cxxopts.hpp"

#include "instance.h"
#include "bnb_mapfts_cplex.h"
#include "bpc_tree.h"
#include "cbs_ts.h"
#include "bnb_lag.h"

void solve_lp_cplex(shared_ptr<Instance> instance);
void solve_mip_cplex(shared_ptr<Instance> instance);
void solve_lp_colgen(shared_ptr<Instance> instance);
void solve_mip_bpc(shared_ptr<Instance> instance);
void solve_cbs(shared_ptr<Instance> instance);
void solve_lbb(shared_ptr<Instance> instance);

Param parse_args(int argc, char **argv);

int main(int argc, char **argv)
{
  debugln("{}\n", string(60, '*'));
  debugln("In DEBUG mode");

  Param param = parse_args(argc, argv);
  debugln("Done: load parameters");

  auto instance = std::make_shared<Instance>(param);
  debugln("Done: load the instance");
#if (MSG_LEVEL >= 1)
  instance->sylog.print_inst_summary();
#endif

  debugln("{}", string(56, '*'));
  debugln("Algorithm Logs");
  debugln("{}\n", string(56, '~'));

  switch (param.algorithm)
  {
  case Param::Alg::lp_cplex:
    solve_lp_cplex(instance);
    break;
  case Param::Alg::mip_cplex:
    solve_mip_cplex(instance);
    break;
  case Param::Alg::ccg:
    solve_lp_colgen(instance);
    break;
  case Param::Alg::bpc:
    solve_mip_bpc(instance);
    break;
  case Param::Alg::cbs:
    solve_cbs(instance);
    break;
  case Param::Alg::lbb:
    solve_lbb(instance);
    break;
  default:
    break;
  }
  debugln("{}", string(56, '~'));
  debugln("Ends normally");
  debugln("{}\n", string(56, '*'));
  return 0;
}

Param parse_args(int argc, char **argv)
{
  Param param;

  // default values
  string default_file       = "../data/movingai_2019/brc202d-random-22.scen";
  string default_num_agent  = "30";
  string default_formul     = "R";
  string default_ctrs       = "1234";
  string default_time_limit = "600";
  string default_gap_limit  = "0";
  string default_ts         = "2";
  string default_alg        = "2"; // 0:lp  1:mip  2:ccg  3:bcp  4:cbs  5:lbb
  string default_subalg     = "1"; // 0:Dijkstar  1:A*
  string default_stepsize   = "6"; // 0: 0.001  1: 0.005  2: 0.01  3: 0.05  4: 0.1  5: 1/iter_count  6: polyak
  string default_cut_rounds = "5";
  string default_cut_depth  = "-1";

  /* Parsing shell parameters */
  // creating program options using cxxopts
  cxxopts::Options cmdline_options(argv[0],
                                   "MAPF-TS - Solution approches for the MAPFP with time-spacing constraints");
  // cmdline_options.show_positional_help();
  cmdline_options.positional_help("instance_file").show_positional_help();
  cmdline_options.add_options()
  ("help", "Print help")
  ("f,file", "Path to instance file", cxxopts::value<vector<string>>()->default_value(default_file))
  ("n,num_agent", "Read first N agents only", cxxopts::value<int>()->default_value(default_num_agent))
  ("o,formulation", "Read formulation options :\n\tV: vertex\n\tA: arc\n\tW: wait_vertex\n\tR: time_range_outgoing", cxxopts::value<string>()->default_value(default_formul))
  ("c,constraints", "Read constraint options :\n\tV: vertex\n\tA: arc\n\tW: wait_vertex\n\tR: time_range_outgoing\n\tO: old_time_spacing\n\tT: time_spacing\n\tC: clique\n\t1: TVC1\n\t2: TVC2\n\t3: TVC3\n\t4: TVC4", cxxopts::value<string>()->default_value(default_ctrs))
  ("t,time-limit", "Time limit in seconds", cxxopts::value<double>()->default_value(default_time_limit))("g,gap-limit", "Solve to an optimality gap", cxxopts::value<double>()->default_value(default_gap_limit))
  ("k,time-spacing", "Minimum time-spacing param", cxxopts::value<int>()->default_value(default_ts))
  ("a,algorithm", "Algorithm to use :\n\t0: lp_cplex\n\t1: mip_cplex\n\t2: ccg\n\t3: bcp\n\t4: cbs\n\t5: lbb", cxxopts::value<int>()->default_value(default_alg))
  ("s,sub-algorithm", "Sub algorithm for pricing :\n\t0: dijkstra\n\t1: A*", cxxopts::value<int>()->default_value(default_subalg))
  ("r,stepsize-rule", "Step size rule for subgradient algorithm: \n\t0: fixed 0.001\n\t1: fixed 0.005\n\t2: fixed 0.01\n\t3: fixed 0.05\n\t4: fixed 0.1\n\t5: 1/iter_count\n\t6: Polyak", cxxopts::value<int>()->default_value(default_stepsize))
  ("u,cut-rounds", "Number of rounds to separate cuts in each node", cxxopts::value<int>()->default_value(default_cut_rounds))
  ("d,cut-depth", "depth of BB tree to separate cuts", cxxopts::value<int>()->default_value(default_cut_depth));
  cmdline_options.parse_positional({"file"});

  // Parse options.
  auto result = cmdline_options.parse(argc, argv);

  // Print help.
  if (result.count("help"))
  {
    println("{}", cmdline_options.help());
    exit(0);
  }

  // Get path to instance.
  param.instance_file = result["file"].as<vector<string>>().at(0);

  // Get agents limit.
  param.num_agents = result["num_agent"].as<int>();

  // Get formulation
  param.formul_string = result["formulation"].as<string>();

  // Get constraints
  param.ctrs_string = result["constraints"].as<string>();

  // Get time limit.
  param.time_limit = result["time-limit"].as<double>();

  // Get optimality gap limit.
  param.gap_limit = result["gap-limit"].as<double>();

  // Get time-spacing parameter.
  param.time_spacing = result["time-spacing"].as<int>();

  // Get algorithm parameter.
  param.algorithm = (Param::Alg)result["algorithm"].as<int>();

  // Get subalgorithm parameter.
  param.sub_algorithm = (Param::SubAlg)result["sub-algorithm"].as<int>();

  // Get stepsize rule parameter
  param.stepsize_rule = (Param::StepSizeRule)result["stepsize-rule"].as<int>();

  // Get number of cuts to generate in each separation
  param.cut_rounds = result["cut-rounds"].as<int>();

  // Get depth of BB tree to separate cuts
  param.cut_depth = result["cut-depth"].as<int>();

  param.alg_str = Param::strings_alg[param.algorithm];
  param.subalg_str = Param::strings_subalg[param.sub_algorithm];
  param.stepsize_str = Param::strings_stepsize[param.stepsize_rule];

  return param;
}

void solve_lp_cplex(shared_ptr<Instance> instance)
{
  BnBSolver bnb_solver(instance);
  bool solve_linear = true;
  bool show_solution = true;
  bnb_solver.solve(solve_linear, show_solution);
  debugln("[DEBUG] Done: solving LPR via CPLEX");
}

void solve_mip_cplex(shared_ptr<Instance> instance)
{
  BnBSolver bnb_solver(instance);
  bool solve_linear = false;
  bool show_solution = true;
  bnb_solver.solve(solve_linear, show_solution);
  debugln("Done: solving MIP via CPLEX BnB");
}

void solve_lp_colgen(shared_ptr<Instance> instance)
{
  BPCTree bpc_tree(instance);
  debugln("[DEBUG] Done: construct BPCTree");
  bpc_tree.solve_root();
  debugln("[DEBUG] Done: solve BPCTree root");
}

void solve_mip_bpc(shared_ptr<Instance> instance)
{
  BPCTree bpc_tree(instance);
  debugln("[DEBUG] Done: construct BPC Tree");
  bpc_tree.explore_tree();
  debugln("[DEBUG] Done: explore BPC Tree");
}

void solve_cbs(shared_ptr<Instance> instance)
{
  cbs_ts(instance);
  debugln("[DEBUG] Done: running CBS");
}

void solve_lbb(shared_ptr<Instance> instance)
{
  bnb_lag(instance);
  debugln("[DEBUG] Done: running LBB");
}