// #ifndef _SY_BRANCHING_RULE_H_
// #define _SY_BRANCHING_RULE_H_

#include "column.h"

class BranchingRule
{
public:
  virtual void add_erased_arcs(ErasedArcsOfAgents& erased_arcs) = 0;
  virtual bool is_column_compatible(shared_ptr<Column> col) = 0;
  virtual void print() = 0;

  virtual ~BranchingRule() = default;
};

class ForbidArcs : public BranchingRule
{
public:
  int agent_id;
  sy_map<tuple<shared_ptr<Arc>, int>, bool> forbided;

  ForbidArcs(int agent_id, sy_map<tuple<shared_ptr<Arc>, int>, bool> forbided) : agent_id(agent_id), forbided(forbided) {}

  void add_erased_arcs(ErasedArcsOfAgents& erased_arcs);
  bool is_column_compatible(shared_ptr<Column> col);
  void print();
};

// #endif