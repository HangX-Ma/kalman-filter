cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

project(kf LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    message(STATUS "Platform: Linux ")
    set(CMAKE_CXX_COMPILER "/usr/bin/clang++-14")
elseif (CMAKE_SYSTEM_NAME MATCHES "Windows")
    message(STATUS "Platform: Windows")
    set(CMAKE_CXX_COMPILER "E:/mingw64/bin/g++.exe")
else ()
    message(FATAL_ERROR "Platform: Unknown")
endif ()

# generate a compile commands file as complete database for vim-YouCompleteMe or some other similar tools
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
add_custom_target(
    copy-compile-commands ALL
    ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_BINARY_DIR}/compile_commands.json
    ${CMAKE_CURRENT_LIST_DIR}
)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "Please select another folder to build the project.")
endif()

# configurations for vscode CMake tools
if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
    message(STATUS "Maximum optimization for speed")
    add_compile_options(-Ofast)
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
    message(STATUS "Maximum optimization for speed, debug info included")
    add_compile_options(-Ofast -g)
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "MinSizeRel")
    message(STATUS "Maximum optimization for size")
    add_compile_options(-Os)
else()
    message(STATUS "Minimal optimization, debug info included")
    add_compile_options(-Og -g)
endif()

add_compile_options(
    -Wall # print warning info
    -Wextra # print extra warning info
    # -Wconversion # data conversion truncation warning
    # -Wsign-conversion # sign <-> unsigned conversion warning
)

add_definitions(
    "-DKF_USE_MATPLOTLIB"
)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Library target
add_library(matplotlib_cpp INTERFACE)
target_include_directories(matplotlib_cpp
  INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/examples>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(matplotlib_cpp INTERFACE
  cxx_std_17
)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

include_directories(${Python3_INCLUDE_DIR})
target_link_libraries(matplotlib_cpp INTERFACE
  Python3::Python
  Python3::Module
)

find_package(Python3 COMPONENTS NumPy)
if(Python3_NumPy_FOUND)
  target_link_libraries(matplotlib_cpp INTERFACE
    Python3::NumPy
  )
else()
  target_compile_definitions(matplotlib_cpp INTERFACE WITHOUT_NUMPY)
endif()

# KF const mode target
add_executable(kf_const_model ${CMAKE_SOURCE_DIR}/demo/kf_const_model.cc kalman_filter.h)
target_link_libraries(kf_const_model  matplotlib_cpp)
target_include_directories(kf_const_model
    PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/library/matplotlib-cpp
        ${CMAKE_SOURCE_DIR}/library/eigen-3.4.0
)

# KF acceleration noise mode target
add_executable(kf_acc_noise_model ${CMAKE_SOURCE_DIR}/demo/kf_acc_noise_model.cc kalman_filter.h)
target_link_libraries(kf_acc_noise_model  matplotlib_cpp)
target_include_directories(kf_acc_noise_model
    PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/library/matplotlib-cpp
        ${CMAKE_SOURCE_DIR}/library/eigen-3.4.0
)

# KF Visual Servo
add_executable(kf_vs_const_model ${CMAKE_SOURCE_DIR}/demo/kf_vs_const_model.cc kalman_filter.h)
target_link_libraries(kf_vs_const_model matplotlib_cpp)
target_include_directories(kf_vs_const_model
    PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/demo
        ${CMAKE_SOURCE_DIR}/library/matplotlib-cpp
        ${CMAKE_SOURCE_DIR}/library/eigen-3.4.0
)