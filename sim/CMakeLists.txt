add_executable(sim sim.cpp)
add_executable(sim_q sim_q.cpp)
add_executable(sim_s sim_s.cpp)
add_executable(sim_sweep_initial sim_sweep_initial.cpp)
add_executable(sim_sweep_initial_fixed sim_sweep_initial_fixed.cpp)
add_executable(sim_rr sim_rr.cpp)

set(xtensor_DIR "${CMAKE_SOURCE_DIR}/lib64/cmake/xtensor/")
set(xsimd_DIR "${CMAKE_SOURCE_DIR}/lib64/cmake/xsimd/")
set(xtl_DIR "${CMAKE_SOURCE_DIR}/share/cmake/xtl/")

list(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/lib64")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -march=native")

find_package(OpenMP)
find_package(xtensor REQUIRED)
message(STATUS "xtensor_DIR: ${xtensor_DIR}")

set(XTENSOR_USE_XSIMD 1)

target_compile_features(gbp PUBLIC cxx_std_17)

#if(OpenMP_CXX_FOUND)
#    target_link_libraries(sim PUBLIC OpenMP::OpenMP_CXX)
#    target_link_libraries(sim_q PUBLIC OpenMP::OpenMP_CXX)
#    target_link_libraries(sim_s PUBLIC OpenMP::OpenMP_CXX)
#    target_link_libraries(sim_sweep_initial PUBLIC OpenMP::OpenMP_CXX)
#    target_link_libraries(sim_sweep_initial_fixed PUBLIC OpenMP::OpenMP_CXX)
#    target_link_libraries(sim_rr PRIVATE gbp m stdc++ stdc++fs xtensor  OpenMP::OpenMP_CXX)
#    target_link_libraries(sim_sn PRIVATE gbp m stdc++ stdc++fs xtensor  OpenMP::OpenMP_CXX)
#endif()

target_link_libraries(sim PRIVATE gbp m stdc++ stdc++fs xtensor pthread)
target_link_libraries(sim_q PRIVATE gbp m stdc++ stdc++fs xtensor pthread)
target_link_libraries(sim_s PRIVATE gbp m stdc++ stdc++fs xtensor pthread)
target_link_libraries(sim_sweep_initial PRIVATE gbp m stdc++ stdc++fs xtensor pthread)
target_link_libraries(sim_sweep_initial_fixed PRIVATE gbp m stdc++ stdc++fs xtensor pthread)
target_link_libraries(sim_rr PRIVATE gbp m stdc++ stdc++fs xtensor pthread)