cmake_minimum_required(VERSION 3.30)
project(PassiveMagneticAttitudeOrbitSimulator)

# list(APPEND CMAKE_MODULE_PATH "/usr/share/cmake/geographiclib")

find_package(Eigen3 CONFIG REQUIRED)
find_package(GeographicLib CONFIG REQUIRED)
find_package(tomlplusplus CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED COMPONENTS
    program_options
    filesystem
)

add_subdirectory(third-party)

#
# Library
#

add_library(pmaos_core STATIC)

set_target_properties(pmaos_core PROPERTIES
    CXX_SCAN_FOR_MODULES OFF
    INTERPROCEDURAL_OPTIMIZATION_RELEASE ON
)

target_compile_definitions(pmaos_core PRIVATE
    WEATHER_DATA_PATH="${SPACE_WEATHER_FILE_PATH}"
)

target_compile_features(pmaos_core PUBLIC
    cxx_std_23
)

target_compile_options(pmaos_core PUBLIC
    -Wall
    -Wextra
)

target_include_directories(pmaos_core PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/source"
)

target_link_libraries(pmaos_core
    PRIVATE
        Boost::boost # for odeint
        GeographicLib::GeographicLib
        nrlmsise-00::nrlmsise-00
        tomlplusplus::tomlplusplus
    PUBLIC
        Eigen3::Eigen
)

target_sources(pmaos_core PRIVATE
    "source/aos/cli.cpp"
    "source/aos/cli.hpp"
    "source/aos/components/hysteresis_rod.cpp"
    "source/aos/components/hysteresis_rod.hpp"
    "source/aos/components/hysteresis_rods.cpp"
    "source/aos/components/hysteresis_rods.hpp"
    "source/aos/components/permanent_magnet.cpp"
    "source/aos/components/permanent_magnet.hpp"
    "source/aos/components/spacecraft_face.cpp"
    "source/aos/components/spacecraft_face.hpp"
    "source/aos/components/spacecraft_faces.cpp"
    "source/aos/components/spacecraft_faces.hpp"
    "source/aos/components/spacecraft_shape.cpp"
    "source/aos/components/spacecraft_shape.hpp"
    "source/aos/components/spacecraft.cpp"
    "source/aos/components/spacecraft.hpp"
    "source/aos/core/constants.hpp"
    "source/aos/core/state.cpp"
    "source/aos/core/state.hpp"
    "source/aos/core/types.hpp"
    "source/aos/environment/details/environment_impl.cpp"
    "source/aos/environment/details/environment_impl.hpp"
    "source/aos/environment/environment.cpp"
    "source/aos/environment/environment.hpp"
    "source/aos/environment/nrlmsise.cpp"
    "source/aos/environment/nrlmsise.hpp"
    "source/aos/environment/orbital_mechanics.cpp"
    "source/aos/environment/orbital_mechanics.hpp"
    "source/aos/environment/space_weather.cpp"
    "source/aos/environment/space_weather.hpp"
    "source/aos/simulation/config.cpp"
    "source/aos/simulation/config.hpp"
    "source/aos/simulation/dynamics.cpp"
    "source/aos/simulation/dynamics.hpp"
    "source/aos/simulation/observer.cpp"
    "source/aos/simulation/observer.hpp"
    "source/aos/simulation/simulation.cpp"
    "source/aos/simulation/simulation.hpp"
    "source/aos/simulation/spacecraft_dynamics.cpp"
    "source/aos/simulation/spacecraft_dynamics.hpp"
    "source/aos/verify/hysteresis_loop_dynamics.cpp"
    "source/aos/verify/hysteresis_loop_dynamics.hpp"
    "source/aos/verify/hysteresis_observer.cpp"
    "source/aos/verify/hysteresis_observer.hpp"
    "source/aos/verify/hysteresis.cpp"
    "source/aos/verify/hysteresis.hpp"
    "source/aos/verify/verification_observer.cpp"
    "source/aos/verify/verification_observer.hpp"
)

#
# Executables
#

add_executable(pmaos_run "source/main.cpp")
set_target_properties(pmaos_run PROPERTIES CXX_SCAN_FOR_MODULES OFF INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
target_link_libraries(pmaos_run PRIVATE pmaos_core)

add_executable(pmaos_vh "source/verify_hysteresis.cpp")
set_target_properties(pmaos_vh PROPERTIES CXX_SCAN_FOR_MODULES OFF INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
target_link_libraries(pmaos_vh PRIVATE pmaos_core)

add_executable(pmaos_vs "source/verify_simulation.cpp")
set_target_properties(pmaos_vs PROPERTIES CXX_SCAN_FOR_MODULES OFF INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
target_link_libraries(pmaos_vs PRIVATE pmaos_core)
