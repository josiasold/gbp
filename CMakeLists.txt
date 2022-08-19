cmake_minimum_required(VERSION 3.16) 

project(gbp
        VERSION 1.0
        LANGUAGES CXX)


# compiled library code
add_subdirectory(src)

# simulation executables
add_subdirectory(sim)


# "installation"
install(TARGETS sim
        DESTINATION ${PROJECT_SOURCE_DIR}/sim)
install(TARGETS sim_q
        DESTINATION ${PROJECT_SOURCE_DIR}/sim)
install(TARGETS sim_s
        DESTINATION ${PROJECT_SOURCE_DIR}/sim)
install(TARGETS sim_sweep_initial
        DESTINATION ${PROJECT_SOURCE_DIR}/sim)
install(TARGETS sim_sweep_initial_fixed
        DESTINATION ${PROJECT_SOURCE_DIR}/sim)
install(TARGETS sim_rr
        DESTINATION ${PROJECT_SOURCE_DIR}/sim)