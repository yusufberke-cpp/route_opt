cmake_minimum_required(VERSION 3.20)

project(sr LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)

add_executable(sr
    main.cpp
    dp_opt.cpp
    route_organizer.cpp
)

target_include_directories(sr PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    C:/vcpkg/installed/x64-windows/include
)

target_link_libraries(sr PRIVATE
    OpenMP::OpenMP_CXX
)

target_include_directories(sr PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/external
)
