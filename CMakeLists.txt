cmake_minimum_required(VERSION 2.6)

#######################
# project information #
#######################

project(bmrm)
set (BMRM_VERSION_MAJOR 2)
set (BMRM_VERSION_MINOR 1)

#####################
# required packages #
#####################

find_package(Boost 1.42)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
    message(STATUS "Boost found.")
else()
    message(STATUS "Boost *NOT* found.")
endif()

# define core libraries
set(CORELIBS ${Boost_LIBRARIES})

set(CPLEX_BASE_PATH "~/local" CACHE PATH "Additional install path for cplex")
set(CONCERT_BASE_PATH "~/local" CACHE PATH "Additional install path for concert")
set(CPLEX_SYSTEM "x86-64_sles10_4.1" CACHE PATH "Specify the system of cplex")
set(CPLEX_LIBFORMAT "static_pic" CACHE PATH "Specify the library format of cplex")

add_definitions(-DHAVE_CPLEX)
add_definitions(-DIL_STD)

include_directories(${CPLEX_BASE_PATH}/include)
include_directories(${CONCERT_BASE_PATH}/include)

link_directories(${CPLEX_BASE_PATH}/lib/${CPLEX_SYSTEM}/${CPLEX_LIBFORMAT})
link_directories(${CONCERT_BASE_PATH}/lib/${CPLEX_SYSTEM}/${CPLEX_LIBFORMAT})

#######################
# project directories #
#######################

include_directories(${PROJECT_SOURCE_DIR}/data)
include_directories(${PROJECT_SOURCE_DIR}/linesearch)
include_directories(${PROJECT_SOURCE_DIR}/loss)
include_directories(${PROJECT_SOURCE_DIR}/model)
include_directories(${PROJECT_SOURCE_DIR}/regularizer)
include_directories(${PROJECT_SOURCE_DIR}/solver)
include_directories(${PROJECT_SOURCE_DIR}/solver/bmrminnersolver)
include_directories(${PROJECT_SOURCE_DIR}/tools)
include_directories(${PROJECT_SOURCE_DIR}/utilities)

################
# source files #
################

file (GLOB src_data ${PROJECT_SOURCE_DIR}/data/*.cpp)
file (GLOB src_linesearch ${PROJECT_SOURCE_DIR}/linesearch/*.cpp)
file (GLOB src_loss ${PROJECT_SOURCE_DIR}/loss/*.cpp)
file (GLOB src_model ${PROJECT_SOURCE_DIR}/model/*.cpp)
file (GLOB src_regularizer ${PROJECT_SOURCE_DIR}/regularizer/*.cpp)
file (GLOB src_tools ${PROJECT_SOURCE_DIR}/tools/*.cpp)
file (GLOB src_utilities ${PROJECT_SOURCE_DIR}/utilities/*.cpp)
file (GLOB src_solver ${PROJECT_SOURCE_DIR}/solver/bmrminnersolver/l2n2*.cpp)
set (src_solver ${src_solver};${PROJECT_SOURCE_DIR}/solver/bmrminnersolver/bmrminnersolver.cpp)
set (src_solver ${src_solver};${PROJECT_SOURCE_DIR}/solver/bmrm.cpp)
if (BT_FOUND)
    set (src_solver ${src_solver};${PROJECT_SOURCE_DIR}/solver/bt.cpp)
endif()
if (Clp_FOUND)
    set (src_solver ${src_solver};${PROJECT_SOURCE_DIR}/solver/bmrminnersolver/l1n1_clp.cpp)
endif()

#####################
# compiler settings #
#####################

set (CMAKE_CXX_FLAGS_RELEASE "-O3 -fomit-frame-pointer")
set (CMAKE_CXX_FLAGS_DEBUG   "-g")

############
# binaries #
############

add_subdirectory(linear-bmrm)
