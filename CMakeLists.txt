cmake_minimum_required(VERSION 3.30)
project(AttitudeOrbitSimulator)

find_package(Eigen3 CONFIG REQUIRED)
find_package(GeographicLib CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED COMPONENTS
    program_options
)

add_executable(AttitudeOrbitSimulator)

target_compile_features(AttitudeOrbitSimulator PRIVATE
    cxx_std_23
)

target_compile_options(AttitudeOrbitSimulator PRIVATE
    -Wall
    -Wextra
)

target_include_directories(AttitudeOrbitSimulator PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(AttitudeOrbitSimulator PRIVATE
    Eigen3::Eigen
    Boost::program_options
    GeographicLib::GeographicLib
)

target_sources(AttitudeOrbitSimulator PRIVATE
    main.cpp
)
