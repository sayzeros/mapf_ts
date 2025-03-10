#ifndef SY_CLIQUE_H
#define SY_CLIQUE_H

#include "graph.h"
#include "column.h"

extern "C"
{
#include "cliquer-1.21/cliquer.h"
}

struct Clique
{
  vector<tuple<int, shared_ptr<Vertex>, shared_ptr<Vertex>, int, shared_ptr<Column>>> iat_list;
};

boolean clique_print_time_with_limit(int level, int i, int n, int max,
			  double cputime, double realtime,
			  clique_options *opts); 

#endif