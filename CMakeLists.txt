cmake_minimum_required(VERSION 3.30)
project(AttitudeOrbitSimulator)

find_package(Eigen3 CONFIG REQUIRED)
find_package(GeographicLib CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED COMPONENTS
    program_options
)

add_executable(AttitudeOrbitSimulator)

set_target_properties(AttitudeOrbitSimulator PROPERTIES
    CXX_SCAN_FOR_MODULES OFF
)

target_compile_features(AttitudeOrbitSimulator PRIVATE
    cxx_std_23
)

target_compile_options(AttitudeOrbitSimulator PRIVATE
    -Wall
    -Wextra
)

target_include_directories(AttitudeOrbitSimulator PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/source"
)

target_link_libraries(AttitudeOrbitSimulator PRIVATE
    Eigen3::Eigen
    Boost::program_options
    GeographicLib::GeographicLib
)

target_sources(AttitudeOrbitSimulator PRIVATE
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
    "source/aos/simulation/dynamics.cpp"
    "source/aos/simulation/dynamics.hpp"
    "source/aos/verify/hysteresis.cpp"
    "source/aos/verify/hysteresis.hpp"
    "source/main.cpp"
)
