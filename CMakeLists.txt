cmake_minimum_required(VERSION 3.16)
project(LargeKernelMatting VERSION 0.8
                           LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(JPEG REQUIRED)
find_package(Eigen3 REQUIRED)

file(GLOB SOURCES
        src/io/*.cpp
        src/matting/*.cpp
        )

file(GLOB HEADERS
        include/io/*.h
        include/matting/*.h
        include/timer/*.h
        )

add_executable(matting main.cpp ${SOURCES} ${HEADERS})
target_link_libraries(matting JPEG::JPEG Eigen3::Eigen)
target_compile_options(matting PRIVATE -Wall -Wextra)
