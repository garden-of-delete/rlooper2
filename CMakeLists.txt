cmake_minimum_required(VERSION 3.6)
project(rlooper2)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES structure.h structure.cpp biophysics.h gene.h Rloop_equilibrium_model.h gene.cpp Rloop_equilibrium_model.cpp model.h model.cpp exception_handling.h simulation.h simulation.cpp windower.cpp windower.h lumberjack.cpp lumberjack.h Rloop_dynamic_model.cpp Rloop_dynamic_model.h)
set(ENSEMBLE_ANALYZER_SOURCE_FILES ensemble_analyzer_main.cpp)
set(PERLOOP_SOURCE_FILES perloop_main.cpp)
set(EXLECTED_LENGTH_SOURCE_FILES expected_length_main.cpp)
set(DYNAMIC_SIM_SOURCE_FILES dynamic_sim_main.cpp)

add_executable(ensemble_analyzer ${SOURCE_FILES} ${ENSEMBLE_ANALYZER_SOURCE_FILES})
add_executable(perloop ${SOURCE_FILES} ${PERLOOP_SOURCE_FILES})
add_executable(dynamic_sim ${SOURCE_FILES} ${DYNAMIC_SIM_SOURCE_FILES})