file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${gbp_SOURCE_DIR}/include/gbp/*.hpp")

message(STATUS "bp_gbp_SOURCE_DIR: ${bp_gbp_SOURCE_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
message(STATUS "HEADER_LIST: ${HEADER_LIST}")

add_library(gbp error_channels.cpp timing.cpp factor.cpp exceptions.cpp util.cpp varset.cpp region_graph.cpp decoder.cpp ${HEADER_LIST})

target_include_directories(gbp PUBLIC ../include)
target_include_directories(gbp PUBLIC ../lib)
target_include_directories(gbp PUBLIC ../lib64)

target_include_directories(gbp PUBLIC /usr/local/include)

find_library(
    GMP
    NAMES gmp libgmp
    HINTS "/usr/local/lib"
)

if (NOT GMP)
    message(SEND_ERROR "Did not find lib gmp")
endif ()
message(STATUS "GMP: ${GMP_INCLUDE_DIRS}")

find_library(
    NTL
    NAMES ntl libntl
    HINTS "${gbp_SOURCE_DIR}/lib"
)

# if (NOT NTL)
#     message(SEND_ERROR "Did not find lib ntl")
# endif ()

find_library(
    LEMON
    NAMES lemon liblemon emon libemon
    HINTS "${CMAKE_SOURCE_DIR}/lib"
)

if (NOT LEMON)
    message(SEND_ERROR "Did not find lib lemon")
endif ()


set(XTENSOR_USE_XSIMD 1)

list(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/lib64")

set(xtensor_DIR "${CMAKE_SOURCE_DIR}/lib64/cmake/xtensor/")
set(xsimd_DIR "${CMAKE_SOURCE_DIR}/lib64/cmake/xsimd/")
set(xtl_DIR "${CMAKE_SOURCE_DIR}/share/cmake/xtl/")

message(STATUS "xtensor_DIR: ${xtensor_DIR}")

find_package(
    xtensor REQUIRED
)

target_compile_features(gbp PUBLIC cxx_std_17)


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -march=native")
# SET(PLATFORM_SPECIFIC_LIBS "-lpthread")


target_link_libraries(gbp PUBLIC ${LEMON} stdc++ stdc++fs ${NTL} ${GMP} m xtensor)

source_group(
    TREE "${PROJECT_SOURCE_DIR}/include"
    PREFIX "Header Files" 
    FILES ${HEADER_LIST}
)