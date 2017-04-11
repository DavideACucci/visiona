cmake_minimum_required(VERSION 2.8)

project(detect)

## dependencies

include_directories(/usr/include/eigen3)

include_directories(include)

## detectlib

set(CMAKE_CXX_FLAGS "-std=c++11 -g -o3 ${OTHERS_DEFINE}")

add_library(visiona SHARED src/Visiona.cpp src/DebugPlotConfig.cpp src/MarkerDetector_impl.cpp src/MarkerDetectorConfig.cpp src/Timer.cpp)

target_link_libraries(visiona config++ opencv_core opencv_calib3d opencv_highgui opencv_imgproc)

## batchdetector

add_executable(batchdetector src/batchdetector.cpp)

target_link_libraries(batchdetector visiona)

## install

install(TARGETS visiona DESTINATION lib)

install(FILES include/Visiona.h DESTINATION include)
install(FILES include/DebugPlotConfig.h DESTINATION include)
install(FILES include/MarkerDetectorConfig.h DESTINATION include)