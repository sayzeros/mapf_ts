#include <tuple>
#include "tuple_hash.h"
#include "robin_hood.h"

typedef std::tuple<std::string, int, int> MyTuple;

struct hashMyTuple : public std::unary_function<MyTuple, std::size_t>
{
 std::size_t operator()(const MyTuple& k) const
 {
   return hashValue(k);
 }
};

  robin_hood::unordered_map<MyTuple, double, hashMyTuple> faster_map;