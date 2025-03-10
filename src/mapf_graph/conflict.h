#ifndef SY_CONFLICT_H
#define SY_CONFLICT_H

#include "graph.h"

class VertexConflict
{
public:
  int i, j;
  Coord coord;
  int time;

  VertexConflict() : i(-1), j(-1), coord(Coord(-1,-1)), time(-1) {}
  VertexConflict(int i, int j, int x, int y, int time) : i(i), j(j), coord(Coord(x, y)), time(time) {}

  string name() { return fmt::format("({},{},{},{})", i, j, coord.name(), time); }
};

class ArcConflict
{
public:
  int i, j;
  Coord coord_v, coord_w;
  int time;

  ArcConflict() : i(-1), j(-1), coord_v(Coord(-1,-1)), coord_w(Coord(-1,-1)), time(-1) {}
  ArcConflict(int i, int j, int v_x, int v_y, int w_x, int w_y, int time) : i(i), j(j), coord_v(Coord(v_x, v_y)), coord_w(Coord(w_x, w_y)), time(time) {}

  string name() { return fmt::format("({},{},{},{},{})", i, j, coord_v.name(), coord_w.name(), time); }
};

class TimeConflict
{
public:
  int i, j;
  Coord coord;
  int time_i, time_j;

  TimeConflict() : i(-1), j(-1), coord(Coord(-1,-1)), time_i(-1), time_j(-1) {}
  TimeConflict(int i, int j, int x, int y, int time_i, int time_j) : i(i), j(j), coord(Coord(x, y)), time_i(time_i), time_j(time_j) {}

  string name() { return fmt::format("({},{},{},{},{})", i, j, coord.name(), time_i, time_j); }
};

#endif