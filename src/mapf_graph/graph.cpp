#include "graph.h"

int Vertex::counter = 0;
int Arc::counter = 0;

shared_ptr<Vertex> Graph::add_vertex(Coord coord)
{
  auto v = make_shared<Vertex>(coord);
  vertices.push_back(v);
  vertex_by_coord[make_tuple(coord.x, coord.y)] = v;
  vertex_by_id[v->id] = v;
  return v;
}

void Graph::add_arc(shared_ptr<Vertex> from, shared_ptr<Vertex> to)
{
  auto a = make_shared<Arc>(from, to);
  if (from->type == Vertex::Type::NORMAL && to->type == Vertex::Type::NORMAL)
  {
    a->type = Arc::Type::NORMAL;
    a->cost = 1;
  }
  else
  {
    a->type = Arc::Type::DUMMY;
    a->cost = 0;
  }
  arcs.push_back(a);
  if (from != to)
  {
    arcs_outgoing_from[from].push_back(a);
    arcs_incoming_to[to].push_back(a);
  }
  arc_by_from_to[make_tuple(from, to)] = a;
  arc_by_id[a->id] = a;
}

void Graph::add_dummy_vertex(Vertex::Type type, Coord coord, int agent_id)
{
  auto v = make_shared<Vertex>(coord);
  v->type = type;
  v->agent_id = agent_id;
  vertices.push_back(v);
  vertex_by_id[v->id] = v;
  if (type == Vertex::Type::DUMMY_ORIGIN)
  {
    auto real_v = get_vertex(coord);
    add_arc(v, real_v);
    dummy_origs.push_back(make_tuple(v, real_v));
  }
  else
  {
    auto real_v = get_vertex(coord);
    add_arc(real_v, v);
    dummy_dests.push_back(make_tuple(v, real_v));
  }
}

shared_ptr<Vertex> Graph::get_vertex(Coord coord)
{
  return vertex_by_coord[make_tuple(coord.x, coord.y)];
}

shared_ptr<Vertex> Graph::get_vertex(Vertex::Type vtype, int agent_id)
{
  if (vtype == Vertex::Type::DUMMY_ORIGIN)
    return std::get<1>(dummy_origs[agent_id]);
  else
    return std::get<1>(dummy_dests[agent_id]);
}

shared_ptr<Vertex> Graph::get_vertex(int vertex_id)
{
  return vertex_by_id[vertex_id];
}

shared_ptr<Vertex> Graph::get_dummy_vertex(Vertex::Type vtype, int agent_id)
{
  if (vtype == Vertex::Type::DUMMY_ORIGIN)
    return std::get<0>(dummy_origs[agent_id]);
  else
    return std::get<0>(dummy_dests[agent_id]);
}

vector<shared_ptr<Arc>> Graph::get_out_arcs(shared_ptr<Vertex> v)
{
  if (v->type == Vertex::Type::NORMAL)
  {
    vector<shared_ptr<Arc>> trt(arcs_outgoing_from[v]);
    trt.push_back(get_arc(v, v));
    return trt;
  }
  return arcs_outgoing_from[v];
}
vector<shared_ptr<Arc>> Graph::get_in_arcs(shared_ptr<Vertex> v)
{
  if (v->type == Vertex::Type::NORMAL)
  {
    vector<shared_ptr<Arc>> trt(arcs_incoming_to[v]);
    trt.push_back(get_arc(v, v));
    return trt;
  }
  return arcs_incoming_to[v];
}
vector<shared_ptr<Arc>> Graph::get_outout_arcs(shared_ptr<Vertex> v)
{
  return arcs_outgoing_from[v];
}
vector<shared_ptr<Arc>> Graph::get_inin_arcs(shared_ptr<Vertex> v)
{
  return arcs_incoming_to[v];
}

shared_ptr<Arc> Graph::get_arc(shared_ptr<Vertex> from, shared_ptr<Vertex> to)
{
  return arc_by_from_to[make_tuple(from, to)];
}

shared_ptr<Arc> Graph::get_arc(int id)
{
  return arc_by_id[id];
}

string Graph::to_string()
{
  string trt = "";
  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x++)
    {
      if (get_vertex(Coord(x, y)) != nullptr)
        trt += "o ";
      else
        trt += "x ";
    }
    trt += "\n";
  }
  return trt;
}

vector<shared_ptr<Arc>> Graph::get_arc_path(vector<shared_ptr<Vertex>> v_path)
{
  vector<shared_ptr<Arc>> a_path;
  int idx = -1;
  for (auto& v: v_path) {
    idx += 1;
    if (idx >= v_path.size()-1) break;
    a_path.push_back(get_arc(v, v_path[idx+1]));
  }
  return a_path;
}