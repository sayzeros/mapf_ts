#ifndef SY_GRAPH_H
#define SY_GRAPH_H

#include "sy_abbreviate.h"
#include "sy_unordered_map.h"
#include "sy_log.h"

struct Coord
{
  int x;
  int y;

  Coord() {}
  Coord(int x, int y) : x(x), y(y) {}

  int dist_from(Coord other)
  {
    return abs(other.x - x) + abs(other.y - y);
  }

  string name()
  {
    return fmt::format("({},{})", x, y);
  }
};

class Vertex
{
public:
  enum Type
  {
    NORMAL,
    DUMMY_ORIGIN,
    DUMMY_DESTINATION
  };

  int id = -1;
  Type type;
  Coord coord;
  int agent_id = -1; // valid if type is dummy_

  Vertex() {}

  Vertex(Coord coord) : id{++counter}, type(Type::NORMAL), coord(coord) {}

  bool operator<(const Vertex &rhs) const // need to use std::set
  {
    return id < rhs.id;
  }

  string name()
  {
    string trt = "";
    if (type == Type::NORMAL)
      trt += coord.name();
    else if (type == Type::DUMMY_ORIGIN)
    {
      trt += "(o" + to_string(agent_id) + ")";
    }
    else if (type == Type::DUMMY_DESTINATION)
    {
      trt += "(d" + to_string(agent_id) + ")";
    }

    return trt;
  }

private:
  static int counter;
};

class Arc
{
public:
  enum Type
  {
    NORMAL, // normal arcs
    DUMMY   // starting from dummy_origin or end at dummy_destination
  };
  int id = -1;
  shared_ptr<Vertex> from;
  shared_ptr<Vertex> to;
  double cost;
  Type type;

  Arc(shared_ptr<Vertex> from, shared_ptr<Vertex> to) : id{counter++}, from(from), to(to) {}

  string name()
  {
    string trt = "";
    trt += from->name() + "~" + to->name();
    return trt;
  }

private:
  static int counter;
};

using ErasedArcs = sy_map<tuple<shared_ptr<Arc>, int>, bool>;
using ErasedArcsOfAgents = sy_map<int, ErasedArcs>;

/* Augmented directed graph of MAPF(-TS) instance (grid with holes) */
class Graph
{
public:
  int width;
  int height;
  vector<shared_ptr<Vertex>> vertices;
  vector<shared_ptr<Arc>> arcs;

  Graph(){};

  shared_ptr<Vertex> add_vertex(Coord coord);
  void add_arc(shared_ptr<Vertex> from, shared_ptr<Vertex> to);
  void add_dummy_vertex(Vertex::Type type, Coord coord, int agent_id);

  shared_ptr<Vertex> get_vertex(Coord coord);
  shared_ptr<Vertex> get_vertex(Vertex::Type vtype, int agent_id);
  shared_ptr<Vertex> get_vertex(int vertex_id);
  shared_ptr<Vertex> get_dummy_vertex(Vertex::Type vtype, int agent_id);

  vector<shared_ptr<Arc>> get_out_arcs(shared_ptr<Vertex> v);
  vector<shared_ptr<Arc>> get_in_arcs(shared_ptr<Vertex> v);
  vector<shared_ptr<Arc>> get_outout_arcs(shared_ptr<Vertex> v);
  vector<shared_ptr<Arc>> get_inin_arcs(shared_ptr<Vertex> v);

  shared_ptr<Arc> get_arc(shared_ptr<Vertex> from, shared_ptr<Vertex> to);
  shared_ptr<Arc> get_arc(int id);

  string to_string();
  // void print_info();
  vector<shared_ptr<Arc>> get_arc_path(vector<shared_ptr<Vertex>> v_path);

private:
  sy_map<shared_ptr<Vertex>, vector<shared_ptr<Arc>>> arcs_outgoing_from, arcs_incoming_to;

  vector<tuple<shared_ptr<Vertex>, shared_ptr<Vertex>>> dummy_origs, dummy_dests;

  sy_map<tuple<int, int>, shared_ptr<Vertex>> vertex_by_coord;
  sy_map<tuple<int>, shared_ptr<Vertex>> vertex_by_id;
  sy_map<tuple<shared_ptr<Vertex>, shared_ptr<Vertex>>, shared_ptr<Arc>> arc_by_from_to;
  sy_map<int, shared_ptr<Arc>> arc_by_id;
};

#endif