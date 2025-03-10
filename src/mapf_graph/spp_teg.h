#ifndef SY_TEG_H
#define SY_TEG_H

#include "instance.h"
#include "graph.h"
#include "master_solver.h"
#include "bnb_lag.h"

double get_global_altered_cost(int agent_id, shared_ptr<MasterResult> mas_result);

double get_local_altered_cost(int agent_id, shared_ptr<MasterResult> mas_result, shared_ptr<Arc> a, int t);

double solve_spp_dijkstra(Graph &graph, shared_ptr<Vertex> vertex, shared_ptr<Vertex> goal);

tuple<
    vector<shared_ptr<Vertex>>, double, double>
solve_spp_teg_dijkstra(Graph &graph,
                       int agent_id,
                       int num_agent,
                       int ts,
                       shared_ptr<MasterResult> mas_result,
                       ErasedArcs &erased_arcs);

tuple<
    vector<shared_ptr<Vertex>>, double, double>
solve_spp_teg_dijkstra(Graph &graph,
                       int agent_id,
                       int num_agent,
                       int ts,
                       shared_ptr<MasterResult> mas_result,
                       ErasedArcs &erased_arcs,
                       sy_map<tuple<shared_ptr<Vertex>, int>, bool> &vt_occupied);

tuple<
    vector<shared_ptr<Vertex>>, double, double>
solve_spp_teg_astar(Graph &graph,
                    int agent_id,
                    int num_agent,
                    int ts,
                    shared_ptr<MasterResult> mas_result,
                    ErasedArcs &erased_arcs);

tuple<
    vector<shared_ptr<Vertex>>, double, double>
solve_spp_teg_astar(Graph &graph,
                    int agent_id,
                    int num_agent,
                    int ts,
                    shared_ptr<MasterResult> mas_result,
                    ErasedArcs &erased_arcs,
                    sy_map<tuple<shared_ptr<Vertex>, int>, bool> &vt_occupied);

tuple<
    vector<shared_ptr<Vertex>>, double, double>
solve_spp_teg_astar(Graph &graph,
                    int agent_id,
                    int num_agent,
                    int ts,
                    set<constraint_old_lag> ctrs,
                    shared_ptr<MasterResult> mas_result,
                    sy_map<tuple<shared_ptr<Vertex>, int>, bool> &vt_occupied);


#endif