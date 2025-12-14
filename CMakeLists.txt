cmake_minimum_required(VERSION 3.30)
project(AttitudeOrbitSimulator)

find_package(Eigen3 CONFIG REQUIRED)
find_package(GeographicLib CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED COMPONENTS
    program_options
    filesystem
)

#
# Library
#

add_library(attitude_orbit_simulator_core STATIC)

set_target_properties(attitude_orbit_simulator_core PROPERTIES
    CXX_SCAN_FOR_MODULES OFF
)

target_compile_features(attitude_orbit_simulator_core PUBLIC
    cxx_std_23
)

target_compile_options(attitude_orbit_simulator_core PUBLIC
    -Wall
    -Wextra
)

target_include_directories(attitude_orbit_simulator_core PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/source"
)

target_link_libraries(attitude_orbit_simulator_core PRIVATE
    Boost::filesystem
    Boost::program_options
PUBLIC
    Boost::boost
    Eigen3::Eigen
    GeographicLib::GeographicLib
)

target_sources(attitude_orbit_simulator_core PRIVATE
    "source/aos/components/hysteresis_rod.cpp"
    "source/aos/components/hysteresis_rod.hpp"
    "source/aos/components/permanent_magnet.cpp"
    "source/aos/components/permanent_magnet.hpp"
    "source/aos/components/spacecraft.cpp"
    "source/aos/components/spacecraft.hpp"
    "source/aos/core/constants.hpp"
    "source/aos/core/state.hpp"
    "source/aos/core/types.hpp"
    "source/aos/environment/environment.cpp"
    "source/aos/environment/environment.hpp"
    "source/aos/environment/orbital_mechanics.cpp"
    "source/aos/environment/orbital_mechanics.hpp"
    "source/aos/simulation/config.cpp"
    "source/aos/simulation/config.hpp"
    "source/aos/simulation/dynamics.cpp"
    "source/aos/simulation/dynamics.hpp"
    "source/aos/simulation/observer.cpp"
    "source/aos/simulation/observer.hpp"
    "source/aos/simulation/simulation.cpp"
    "source/aos/simulation/simulation.hpp"
    "source/aos/verify/bh_observer.cpp"
    "source/aos/verify/bh_observer.hpp"
    "source/aos/verify/hysteresis_loop_dynamics.cpp"
    "source/aos/verify/hysteresis_loop_dynamics.hpp"
    "source/aos/verify/hysteresis.cpp"
    "source/aos/verify/hysteresis.hpp"
)

#
# Executable
#

add_executable(attitude_orbit_simulator "source/main.cpp")

set_target_properties(attitude_orbit_simulator PROPERTIES
    CXX_SCAN_FOR_MODULES OFF
)

target_link_libraries(attitude_orbit_simulator PRIVATE
    attitude_orbit_simulator_core
)
