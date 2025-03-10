#include "bnb_lag.h"
#include <cassert>


bnb_node* BnB_Tree::generate_root() {
  root->lb = 0;
  root->parent = NULL;
  root->coll = (collision_t) {-1,-1,-1,-1, nullptr};
  root->ctr = (constraint_old_lag) {-1,-1,nullptr};
  root->idx = num_node++;
  root->depth = 0;

  return root;
}

tuple<bnb_node*, bnb_node*> BnB_Tree::get_branches(bnb_node* n) {
  assert(n->coll.i != -1);

  bnb_node* l = new bnb_node();
  l->lb = n->lb;
  l->parent = n;
  l->idx = num_node++;
  l->depth = n->depth+1;
  l->lamb_v = vt_values(n->lamb_v);
  l->lamb_a = at_values(n->lamb_a);
  l->lamb_o = ivt_values(n->lamb_o);
  l->lamb_r = vt_values(n->lamb_r);
  l->coll = (collision_t) {-1,-1,-1,-1, nullptr};
  l->ctr = (constraint_old_lag) {n->coll.i, n->coll.t_i, n->coll.v};
  l->is_leaf = false;

  bnb_node* r = new bnb_node();
  r->lb = n->lb;
  r->parent = n;
  r->idx = num_node++;
  r->depth = n->depth+1;
  r->lamb_v = vt_values(n->lamb_v);
  r->lamb_a = at_values(n->lamb_a);
  r->lamb_o = ivt_values(n->lamb_o);
  r->lamb_r = vt_values(n->lamb_r);
  r->coll = (collision_t) {-1,-1,-1,-1, nullptr};
  r->ctr = (constraint_old_lag) {n->coll.j, n->coll.t_j, n->coll.v};
  r->is_leaf = false;

  return std::tie(l, r);
}

void BnB_Tree::check_out(bnb_node* n) {
  num_visited += 1;
}
