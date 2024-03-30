cmake_minimum_required(VERSION 3.0.0)
project(guidedbilateral VERSION 0.1.0 LANGUAGES CXX CUDA)
enable_language(CUDA)
include(CTest)
enable_testing()

find_package( OpenCV REQUIRED )
find_package( OpenMP )

include_directories( ${OpenCV_INCLUDE_DIRS} )

add_executable(guidedbilateral_cpu cpu_main.cpp)
target_link_libraries( guidedbilateral_cpu ${OpenCV_LIBS} )
target_link_libraries( guidedbilateral_cpu OpenMP::OpenMP_CXX )

add_executable(guidedbilateral_gpu gpu_main.cu)
target_link_libraries( guidedbilateral_gpu ${OpenCV_LIBS} )

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
