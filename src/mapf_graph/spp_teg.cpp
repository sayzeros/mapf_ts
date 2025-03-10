#include "spp_teg.h"
#include "sy_log.h"

#include <queue>

double get_global_altered_cost(int agent_id, shared_ptr<MasterResult> mas_result)
{
  if (mas_result == nullptr) return 0;
  return -mas_result->cardi_duals[agent_id];
}

double get_local_altered_cost(int agent_id, shared_ptr<MasterResult> mas_result, shared_ptr<Arc> a, int t)
{
  if (mas_result == nullptr) return 0;
  double trt = 0;

  if (mas_result->i_at_cost_altered.contains(agent_id) && mas_result->i_at_cost_altered[agent_id].contains(make_tuple(a, t)))
    trt += mas_result->i_at_cost_altered[agent_id][make_tuple(a, t)];

  if (mas_result->at_cost_altered.contains(make_tuple(a, t)))
    trt += mas_result->at_cost_altered[make_tuple(a, t)];

  return trt;
}