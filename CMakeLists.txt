cmake_minimum_required(VERSION 3.10)
project(HyMinSolver)

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -std=c++11")

if(UNIX AND NOT APPLE)
    message(STATUS "OpenMP exist")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -D_GLIBCXX_PARALLEL -DHAS_OMP")
endif()

include_directories(/usr/local/include . ${PROJECT_SOURCE_DIR}/include)
add_executable(HyMin src/main.cpp)

# For the other tools
add_executable(DataCleaner src/tools/data_cleaner.cpp)
set_target_properties(DataCleaner PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")