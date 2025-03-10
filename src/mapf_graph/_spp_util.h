#ifndef _SY_SPP_UTIL_H_
#define _SY_SPP_UTIL_H_


#include <queue>

// Upper bound of path cost;
// needed to reduce effort to initialize dist[v]s
static const double lsp_spp_big_value = 1e6;

// Need to reconstruct multiple paths by increasing reduced_cost order
template <typename T>
struct PathInfo
{
  T last_visited;
  double reduced_cost;

  PathInfo(T &last_visited, double reduced_cost)
      : last_visited(last_visited), reduced_cost(reduced_cost) {}
};
template <typename T>
struct PathInfoCompare
{
  bool operator()(const PathInfo<T> &n1, const PathInfo<T> &n2) const
  {
    return (n1.reduced_cost > n2.reduced_cost);
  }
};
template <typename T>
using ReducedCostQueue = std::priority_queue<PathInfo<T>, vector<PathInfo<T>>, PathInfoCompare<T>>;

#endif