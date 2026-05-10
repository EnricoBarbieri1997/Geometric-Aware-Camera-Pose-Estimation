cmake_minimum_required(VERSION 3.16)
project(ceres_cylinder_calibration CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Homebrew prefix hint (covers both Apple Silicon and Intel)
list(PREPEND CMAKE_PREFIX_PATH "/opt/homebrew" "/usr/local")

find_package(Ceres REQUIRED)
find_package(nlohmann_json REQUIRED)

add_executable(calibrate main.cpp)
target_link_libraries(calibrate PRIVATE Ceres::ceres nlohmann_json::nlohmann_json)
target_compile_options(calibrate PRIVATE -O2)

add_executable(calibrate_sweep sweep.cpp)
target_link_libraries(calibrate_sweep PRIVATE Ceres::ceres nlohmann_json::nlohmann_json)
target_compile_options(calibrate_sweep PRIVATE -O2)
