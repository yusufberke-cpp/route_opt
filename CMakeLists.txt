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
if(OpenMP_CXX_FOUND)
    target_link_libraries(sr PUBLIC OpenMP::OpenMP_CXX)
endif()
target_include_directories(sr PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
