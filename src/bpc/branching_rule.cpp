#include "branching_rule.h"
#include "sy_log.h"

void ForbidArcs::add_erased_arcs(ErasedArcsOfAgents &erased_arcs)
{
  for (auto &[key, is_forbided] : forbided)
  {
    if (is_forbided)
      erased_arcs[agent_id][key] = true;
  }
}

bool ForbidArcs::is_column_compatible(shared_ptr<Column> col)
{
  if (agent_id != col->agent_id)
    return true;
  int t = -2;
  for (auto &arc : col->instance->graph.get_arc_path(col->path))
  {
    t += 1;
    if (forbided[make_tuple(arc, t)])
      return false;
  }
  return true;
}

void ForbidArcs::print()
{
  println("Erased Arcs:");
  println("  agent {}", agent_id);
  for (auto &[key, value] : forbided)
  {
    if (value)
    println("  arc {} at {}", std::get<0>(key)->name(), std::get<1>(key));
  }
}