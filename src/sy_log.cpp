#include "sy_log.h"
#include "parameters.h"

int SYLog::num_lines = 0;

const string SYLog::instance_entities[] = {
    "inst_file",
    "map",
    "type",
    "id",
    "#agents",
    "t_s",
    "formul",
    "ctrs",
    "time_lim",
    "gap_lim",
    "alg",
    "sub_alg",
    "stepsize_rule",
    "#vert",
    "#vert_real",
    "#vert_dumm",
    "#arc",
    "#arc_real",
    "#arc_dumm"};

const string SYLog::algorithm_entities[] = {
    "stat",
    "time",
    "obj",  
    "lb",
    "ub",
    "root_lb",
    "#col",
    "#row",
    "#node",
    "depth",
    "time_mp",
    "time_sp",
    "time_sepa",
    "time_node",
    "time_br",
    "time_heur",
    "#spp",
    "avg_time_sp",
    "#vc",
    "#ac",
    "#tc",
    "#wc",
    "#rc",
    "#oc",
    "#cc",
    "#1c",
    "#2c",
    "#3c",
    "#4c",
    "time_load",
    "time_solve"};

const string SYLog::bnb_entities[] = {
    "stat",
    "time",
    "obj",
    "lb",
    "ub",
    "#col",
    "#row",
    "#node",
    "time_load",
    "time_solve"};

const string SYLog::cg_entities[] = {};

const string SYLog::bpc_entities[] = {};

bool SYLog::is_entity(EntityType type, string entity)
{
  if (type == EntityType::ALG_ENTITY)
  {
    for (auto ent : algorithm_entities)
    {
      if (ent == entity)
      {
        return true;
      }
    }
  }
  else
  {
    for (auto ent : instance_entities)
    {
      if (ent == entity)
      {
        return true;
      }
    }
  }
  return false;
}

void SYLog::Chunk::add_cell(string entity, string data)
{
  line[entity] = data;
}
string SYLog::Chunk::get_cell(string entity)
{
  if (line[entity] == to_string(DOUBLE_NULL))
    return "-";
  if (line[entity] == "")
    return "NA";
  return line[entity];
}

void SYLog::add_inst_info(string entity, int data)
{
  add_inst_info(entity, to_string(data));
}
void SYLog::add_inst_info(string entity, double data)
{
  add_inst_info(entity, to_string(data));
}
void SYLog::add_inst_info(string entity, string data)
{
  debug_assert(is_entity(EntityType::INST_ENTITY, entity));
  inst_chunk.add_cell(entity, data);
}

void SYLog::add_alg_result(string entity, int data, bool new_line)
{
  add_alg_result(entity, to_string(data), new_line);
}
void SYLog::add_alg_result(string entity, long data, bool new_line)
{
  add_alg_result(entity, to_string(data), new_line);
}
void SYLog::add_alg_result(string entity, double data, bool new_line)
{
  add_alg_result(entity, to_string(data), new_line);
}
void SYLog::add_alg_result(string entity, string data, bool new_line)
{
  if (new_line)
  {
    auto new_chunk = Chunk();
    all_chunks.push_back(active_chunk);
    active_chunk = new_chunk;
  }
  debug_assert(is_entity(EntityType::ALG_ENTITY, entity));
  active_chunk.add_cell(entity, data);
}


void SYLog::print_inst_summary()
{
  println("{}", string(60, '='));
  println("Instance Summary");
  println("{}", string(60, '-'));
  for (auto ent : instance_entities)
  {
    println("{:10} : {}", ent, inst_chunk.get_cell(ent));
  }
  println("{}\n", string(60, '='));
}


void SYLog::print_result()
{
  println("{}", string(60, '='));
  println("Results");
  println("{}", string(60, '-'));
  for (auto ent : algorithm_entities)
  {
    if (active_chunk.get_cell(ent) != "NA")
      println("{:11} : {}", ent, active_chunk.get_cell(ent));
  }
  println("{}\n", string(60, '='));
}


void SYLog::print_header(bool is_csv)
{
  string sep = "  ";
  if (is_csv)
  {
    sep = ",";
    for (auto ent : {"map", "type", "id", "#agents", "t_s"})
    {
      print("{}{}", ent, sep);
    }
  }
  if (inst_chunk.get_cell("alg") == "lp_cplex" || inst_chunk.get_cell("alg") == "mip_cplex")
  {
    for (auto ent : bnb_entities)
    {
      print("{}{}", ent, sep);
    }
  }
  println("");
}


void SYLog::print_line(bool is_csv)
{
  string sep = "  ";
  if (is_csv)
  {
    sep = ",";
    for (auto ent : {"map", "type", "id", "#agents", "t_s"})
    {
      print("{}{}", inst_chunk.get_cell(ent), sep);
    }
  }
  if (inst_chunk.get_cell("alg") == "lp_cplex" || inst_chunk.get_cell("alg") == "mip_cplex")
  {
    for (auto ent : bnb_entities)
    {
      print("{}{}", active_chunk.get_cell(ent), sep);
    }
  }
  println("");
}
