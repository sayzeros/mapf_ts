#include "parameters.h"

#include "sy_log.h"

const string Param::strings_alg[] = {
    "lp_cplex",
    "mip_cplex",
    "ccg",
    "bpc",
    "cbs",
    "lbb"};

const string Param::strings_subalg[] = {
    "Dijkstra",
    "A*"};

const string Param::strings_stepsize[] = {
  "fixed1",
  "fixed2",
  "fixed3",
  "fixed4",
  "fixed5",
  "dimini",
  "polyak"
};