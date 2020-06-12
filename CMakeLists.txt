cmake_minimum_required(VERSION 3.16)
project(CFImageMatting)

set(CMAKE_CXX_STANDARD 17)

find_package(JPEG REQUIRED)
find_package(Eigen3 REQUIRED)

add_executable(cf_matting main.cpp src/cfmatting/ImageMatrix.cpp include/cfmatting/ImageMatrix.h src/cfmatting/ImageWindow.cpp include/cfmatting/ImageWindow.h src/cfmatting/SATMatrix.cpp include/cfmatting/SATMatrix.h)
target_link_libraries(cf_matting JPEG::JPEG Eigen3::Eigen)
target_compile_options(cf_matting PRIVATE -Wall -Wextra)
