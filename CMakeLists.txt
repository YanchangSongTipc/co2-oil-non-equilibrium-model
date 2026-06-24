cmake_minimum_required(VERSION 3.14)
project(CO2_NonBalance LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g")

# Header-only autodiff library
add_library(autodiff INTERFACE)
target_include_directories(autodiff INTERFACE ${CMAKE_SOURCE_DIR}/lib)
target_compile_features(autodiff INTERFACE cxx_std_17)

# SAFT-Mie static library
add_library(saft_mie STATIC src/saft_mie.cpp src/saft_mie.h)
target_include_directories(saft_mie PUBLIC ${CMAKE_SOURCE_DIR}/lib)
target_link_libraries(saft_mie PUBLIC autodiff)

# Step 1: Equilibrium solubility LUT
add_executable(calc_weq_lut src/calc_weq_lut.cpp)
target_link_libraries(calc_weq_lut PRIVATE saft_mie)

# Step 2: Mesoscopic transport (tau_ms) LUT
add_executable(calc_taums src/calc_taums.cpp)
target_link_libraries(calc_taums PRIVATE saft_mie)

# Pure CO2 verification test
add_executable(test_saft src/test.cpp)
target_link_libraries(test_saft PRIVATE saft_mie)
