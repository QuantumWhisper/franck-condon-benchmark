# Franck-Condon I-V Simulation -- C++ Implementation
# Build: cmake -B build && cmake --build build
# Run:   ./build/cpp_benchmark quick
#        ./build/cpp_benchmark default

cmake_minimum_required(VERSION 3.16)
project(fc_simulation CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(nlohmann_json REQUIRED)

# OpenMP: Apple Clang needs Homebrew libomp
if(APPLE)
  execute_process(
    COMMAND brew --prefix libomp
    OUTPUT_VARIABLE LIBOMP_PREFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  if(LIBOMP_PREFIX)
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "${LIBOMP_PREFIX}/lib/libomp.dylib")
  endif()
endif()
find_package(OpenMP REQUIRED)

add_executable(cpp_benchmark
    src/laguerre.cpp
    src/fermi_bose.cpp
    src/digamma.cpp
    src/fc_matrix.cpp
    src/regularized.cpp
    src/cotunneling.cpp
    src/rate.cpp
    src/matrix.cpp
    src/solver.cpp
    src/current.cpp
    src/simulate.cpp
    src/json_io.cpp
    src/plotting.cpp
    src/main.cpp
)

target_compile_options(cpp_benchmark PRIVATE
    -O3 -march=native -Wall -Wextra -Wno-unused-parameter
)

target_link_libraries(cpp_benchmark PRIVATE
    Eigen3::Eigen
    nlohmann_json::nlohmann_json
    OpenMP::OpenMP_CXX
)
