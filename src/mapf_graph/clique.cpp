#include "clique.h"

boolean clique_print_time_with_limit(int level, int i, int n, int max,
			  double cputime, double realtime,
			  clique_options *opts) {
	static float prev_time=100;
	static int prev_i=100;
	static int prev_max=100;
	static int prev_level=0;
	FILE *fp=opts->output;
	// int j;

	if (fp==NULL)
		fp=stdout;

	if (ABS(prev_time-realtime)>0.1 || i==n || i<prev_i || max!=prev_max ||
	    level!=prev_level) {
		prev_time=realtime;
		prev_i=i;
		prev_max=max;
		prev_level=level;
	}
  static float time_limit = 300;
  if (realtime > time_limit) return FALSE;
	return TRUE;
}