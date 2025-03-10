#ifndef SY_LOG_H
#define SY_LOG_H

#include "fmt/format.h"
#include <chrono>

#include "sy_abbreviate.h"
#include "sy_unordered_map.h"

  /////////////////////
 // print & println //
/////////////////////
using fmt::print;
#ifdef DEBUG
#define println(format, ...)                \
  do                                        \
  {                                         \
    fmt::print(format "\n", ##__VA_ARGS__); \
    fflush(stdout);                         \
  } while (false)
#define debugln(format, ...) println("[d] " format, ##__VA_ARGS__)
#define debugnewline(num) print("{}", string(num, '\n'))
#define debug(format, ...) fmt::print("[d] " format, ##__VA_ARGS__);
#else
#define println(format, ...)                \
  do                                        \
  {                                         \
    fmt::print(format "\n", ##__VA_ARGS__); \
  } while (false)
#define debugln(format, ...) \
  {                          \
  }
#define debugnewline(num) \
  {                          \
  }
#define debug(format, ...) \
  {                        \
  }
#endif

  //////////////////
 // assert utils //
//////////////////
#ifdef DEBUG
#define err(format, ...) do { \
    fmt::print(stderr, "Error: " format "\n", ##__VA_ARGS__); \
    fmt::print(stderr, "Function: {}\n", __PRETTY_FUNCTION__); \
    fmt::print(stderr, "File: {}\n", __FILE__); \
    fmt::print(stderr, "Line: {}\n", __LINE__); \
    std::abort(); \
} while (false)
#else
#define err(format, ...) do { \
    fmt::print(stderr, "Error: " format "\n", ##__VA_ARGS__); \
    std::abort(); \
} while (false)
#endif

#define release_assert(condition, ...) do { \
    if (!(condition)) err(__VA_ARGS__); \
} while (false)

#ifdef DEBUG
#define debug_assert(condition) release_assert(condition, "{}", #condition)
#else
#define debug_assert(condition) {}
#endif


  ////////////////
 // time stamp //
////////////////
using time_point = std::chrono::high_resolution_clock::time_point;

inline time_point get_now()
{
  return std::chrono::high_resolution_clock::now();
}

inline double get_elapsed(time_point from, time_point to)
{
  return std::chrono::duration_cast<std::chrono::duration<double>>(to - from).count();
}

inline double get_elapsed(time_point from)
{
  return get_elapsed(from, get_now());
}


  ///////////////////////////////////////////
 // class for fomatted consistent logging //
///////////////////////////////////////////
class SYLog
{
public:
  const static string instance_entities[];
  const static string algorithm_entities[];
  const static string bnb_entities[];
  const static string cg_entities[];
  const static string bpc_entities[];

  void add_inst_info(string entity, int data);
  void add_inst_info(string entity, double data);
  void add_inst_info(string entity, string data);
  void add_alg_result(string entity, int data, bool new_line=false);
  void add_alg_result(string entity, long data, bool new_line=false);
  void add_alg_result(string entity, double data, bool new_line=false);
  void add_alg_result(string entity, string data, bool new_line=false);
  void print_inst_summary();
  void print_result();
  void print_header(bool is_csv=false);
  void print_line(bool is_csv=false);

private:
  class Chunk {
    friend class SYLog;
    sy_map<string, string> line;
    void add_cell(string entity, string data);
    string get_cell(string entity);
    Chunk() {line["idx"] = num_lines++;}
  };
  Chunk inst_chunk, active_chunk;
  vector<Chunk> all_chunks;

  enum EntityType
  {
    ALG_ENTITY,
    INST_ENTITY
  };
  bool is_entity(EntityType type, string entity);
  static int num_lines;
};

#endif