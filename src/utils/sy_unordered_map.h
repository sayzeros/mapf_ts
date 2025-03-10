#ifndef SY_UNORDERED_MAP_H
#define SY_UNORDERED_MAP_H

#include "robin_hood.h"
#include "tuple_hash.h"

template<typename Key, typename Value>
using sy_map = robin_hood::unordered_map<Key, Value, tuple_hash::hash<Key> >;

#endif