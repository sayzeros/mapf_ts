cmake_minimum_required(VERSION 3.22)

add_compile_definitions(_MAC_CMAKE_)
add_compile_options("$<$<CONFIG:DEBUG>:-DDEBUG>")

project(
  mapf_ts
  VERSION 0.1
  DESCRIPTION "Algorithms for MAPF-TS via C++"
  LANGUAGES CXX C)

# Set C version.
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED 1)

# Set C++ version.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED 1)

# Set release build.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(CONFIG:DEBUG)
  set(CMAKE_BUILD_TYPE Debug)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
set(CMAKE_CXX_FLAGS "-Wall") 
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fomit-frame-pointer -march=native -flto")

#---------------------------------------------------------------------
# #####################
# # Add source files ##
# #####################
set(BNB_SOURCES
  ./src/bnb_cplex/bnb_mapfts_cplex.cpp
)

set(BPC_SOURCES
  ./src/bpc/bpc_tree.cpp
  ./src/bpc/bpc_node.cpp
  ./src/bpc/column.cpp
  ./src/bpc/master_solver.cpp
  ./src/bpc/sub_solver.cpp
  ./src/bpc/branching_rule.cpp
)

set(CBS_SOURCES
  ./src/cbs/cbs_ts.cpp
)

set(LBB_SOURCES
  ./src/lbb/bnb_lag.cpp
  ./src/lbb/bnb_tree.cpp
)

set(GRAPH_SOURCES
  ./src/mapf_graph/graph.cpp
  ./src/mapf_graph/solution.cpp
  ./src/mapf_graph/spp_teg.cpp
  ./src/mapf_graph/spp_teg_dijkstra.cpp
  ./src/mapf_graph/spp_teg_astar.cpp
  ./src/mapf_graph/clique.cpp
)

set(MAPFTS_SOURCES
  ./src/main.cpp
  ./src/parameters.cpp
  ./src/instance.cpp
  ./src/sy_log.cpp
)

# Crate target
add_executable(${PROJECT_NAME}
  ${MAPFTS_SOURCES}
  ${BNB_SOURCES}
  ${BPC_SOURCES}
  ${CBS_SOURCES}
  ${GRAPH_SOURCES} 
  ${CLIQUER_SOURCES}
  ${LBB_SOURCES}
)

# Macro to collect all subdirectories
MACRO(SUBDIRLIST result curdir)
  file(GLOB children LIST_DIRECTORIES true ${curdir} ${curdir}/*)
  SET(dirlist "")

  FOREACH(child ${children})
    IF(IS_DIRECTORY ${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()

  SET(${result} ${dirlist})
ENDMACRO()

# Include all subdirectory of src
SUBDIRLIST(SUBDIRS src)

FOREACH(subdir ${SUBDIRS})
  target_include_directories(${PROJECT_NAME} PUBLIC ${subdir})
ENDFOREACH()

# Set compile flag -Wall
target_compile_options(${PROJECT_NAME} PUBLIC -Wall)

# ---------------------------------------------------------------------
# ####################
# # Include CXXOPTS ##
# ####################
target_include_directories(${PROJECT_NAME} PUBLIC ./src/utils/cxxopts/include/)

# ---------------------------------------------------------------------
# #######################################
# # Include fmt and robin-hood-hashing ##
# #######################################
add_subdirectory(src/utils/fmt EXCLUDE_FROM_ALL)
target_include_directories(${PROJECT_NAME} PUBLIC .src/utils/fmt/include/ ./src/utils/robin-hood-hashing/src/include/)
target_link_libraries(${PROJECT_NAME} fmt)

# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# ####################
# # Include cliquer ##
# ####################
set(CLIQUER_SOURCES
    ./src/utils/cliquer-1.21/cliquer.c
    ./src/utils/cliquer-1.21/cliquer.h
    ./src/utils/cliquer-1.21/cliquerconf.h
    ./src/utils/cliquer-1.21/graph.c
    ./src/utils/cliquer-1.21/graph.h
    ./src/utils/cliquer-1.21/misc.h
    ./src/utils/cliquer-1.21/reorder.c
    ./src/utils/cliquer-1.21/reorder.h
    ./src/utils/cliquer-1.21/set.h
    )
add_library(cliquer STATIC ${CLIQUER_SOURCES})
target_link_libraries(${PROJECT_NAME} cliquer)
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# ##################
# # Include Cplex ##
# ##################
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/cmake)
find_package(CPLEX)
include_directories(${CPLEX_CONCERT_INCLUDE_DIR} ${CPLEX_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ilocplex) 

# ---------------------------------------------------------------------
# eof