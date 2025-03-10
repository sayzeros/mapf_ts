#ifndef _SY_UTIL_H_
#define _SY_UTIL_H_

#include "sy_log.h"

template<typename FROM, typename TO>
void map_copy_all(TO& to, FROM& from) 
{
  for (auto const& [key, val]: from) {
    to[key] = val;
  }
}

#endif