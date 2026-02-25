add_library(nrlmsise-00 STATIC)
add_library(nrlmsise-00::nrlmsise-00 ALIAS nrlmsise-00)

target_include_directories(nrlmsise-00 SYSTEM PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_sources(nrlmsise-00 PRIVATE
    nrlmsise-00.c
    nrlmsise-00.h
    nrlmsise-00_data.c
)
