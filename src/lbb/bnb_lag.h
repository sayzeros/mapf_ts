#ifndef MY_LBNB_H
#define MY_LBNB_H

#include <vector>
#include <chrono>
#include <queue>

#include "instance.h"
#include "sy_log.h"
#include "mapfts_alg_tools.h"

struct constraint_old_lag{
  int i;
  int t;
  shared_ptr<Vertex> v;
  bool operator<(const constraint_old_lag& r) const
  {
    if (i != r.i) return i < r.i;
    if (t != r.t) return t < r.t;
    return v < r.v;
  }
};

struct bnb_node{
  bnb_node* parent;
  collision_t coll;
  constraint_old_lag ctr;
  double lb;
  double pqueue_index; 
  int idx;
  int depth;
  bool is_leaf = false;
  vt_values lamb_v;
  at_values lamb_a;
  ivt_values lamb_o;
  vt_values lamb_r;
};

struct bnb_node_lag_cmp{
  bool operator()(bnb_node* t, bnb_node* u){
    return t->lb > u->lb;
  }
};

class BnB_Tree : public std::priority_queue<bnb_node*, vector<bnb_node*>, bnb_node_lag_cmp>{
  public:
    bnb_node* root;
    double global_lb;
    double ub = TRIVIAL_BOUND;
    vector<vector<shared_ptr<Vertex>>> ub_routes;
    int num_node;
    int num_visited;
    
    // original methods
    BnB_Tree(): root(new bnb_node()), global_lb(0), ub(TRIVIAL_BOUND), ub_routes(vector<vector<shared_ptr<Vertex>>>()), num_node(0), num_visited(0){};
    bnb_node* generate_root();
    tuple<bnb_node*, bnb_node*> get_branches(bnb_node* n);
    void check_out(bnb_node* n);


  protected:
    inline void update_pqueue_index(bnb_node* node, const int pqueue_index)
    {
        node->pqueue_index = pqueue_index;
    }
};

void bnb_lag(shared_ptr<Instance> instance);

#endif