#include "sy_log.h"

#include "column.h"
#include <iostream>

template <typename Map>
bool map_equal(Map const &lhs, Map const &rhs);

void Column::make_dummy()
{
  obj_coeff = DUMMY_COLUMN_COST;
  for (int id = 0; id < instance->param.num_agents; id++)
    cardi_coeff[id] = 1;
  is_dummy = true;
}

Column::Column(shared_ptr<Instance> instance, int agent_id, vector<shared_ptr<Vertex>> path, double cost)
    : instance(instance), agent_id(agent_id), path(path)
{
  is_dummy = false;

  obj_coeff = cost;

  cardi_coeff[agent_id] = 1;

  int time = 0;
  for (auto v : path)
  {
    if (v->type == Vertex::Type::NORMAL)
    {
      // vertex conflict constraints in formulation
      if (instance->param.formul_string.find('V') != std::string::npos)
        vertex_conflict_coeff[make_tuple(v->id, time)] = 1;
      // arc conflict constraints in formulation
      if (instance->param.formul_string.find('A') != std::string::npos)
      {
        auto w = path[time + 2];
        if (v->id < w->id)
        {
          auto arc = instance->graph.get_arc(v, w);
          arc_conflict_coeff[make_tuple(arc->id, time)] = 1;
        }
        else if (v->id > w->id)
        {
          auto arc = instance->graph.get_arc(w, v);
          arc_conflict_coeff[make_tuple(arc->id, time)] = 1;
        }
      }
      // wait conflict constraints in formulation
      if (instance->param.formul_string.find('W') != std::string::npos)
      {
        auto w = path[time + 2];
        wait_conflict_coeff[make_tuple(v->id, time)] = 1;
        if (v->id != w->id && w->type == Vertex::Type::NORMAL)
        {
          wait_conflict_coeff[make_tuple(w->id, time)] = 1;
        }
      }
      // range conflict constraints in formulation
      if (instance->param.formul_string.find('R') != std::string::npos)
      {
        auto w = path[time + 2];
        if (v->id != w->id)
        {
          for (int tau = 0; tau <= instance->param.time_spacing; tau++)
            range_conflict_coeff[make_tuple(v->id, time - tau)] += 1;
        }
        else
          range_conflict_coeff[make_tuple(v->id, time - instance->param.time_spacing)] += 1;
      }
      // coefficients of separated cuts are set in sub_solver.cpp, since we don't know the separated cuts here
      time += 1;
    }
  }
}

bool Column::operator==(const Column &other) const
{
  return path.size() == other.path.size() &&
         std::equal(path.begin(), path.end(), other.path.begin());
};

template <typename Map>
bool map_equal(Map const &lhs, Map const &rhs)
{
  return (lhs.size() == rhs.size()) && std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

string Column::name()
{
  string trt = fmt::format("|  obj {}  |  ", obj_coeff);
  if (is_dummy)
  {
    trt += "dummy";
  }
  else
  {
    for (auto v : path)
    {
      trt += fmt::format("{} ", v->name());
    }
  }
  trt += " ||  ";

  return trt;
}

std::ostream &operator<<(std::ostream &out, const Column &c)
{
  // auto &instance = c.instance;
  out << "agent " << c.agent_id << "  |  ";
  for (auto &v : c.path)
  {
    out << v->name() << "  ";
  }
  out << "\n";

  return out;
}
