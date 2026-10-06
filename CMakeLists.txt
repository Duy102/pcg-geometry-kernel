cmake_minimum_required(VERSION 3.20)
project(pcg_geometry_kernel LANGUAGES CXX)
find_package(OpenSSL REQUIRED)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
add_library(pcg STATIC
    src/core.cpp
    src/geometry.cpp
    src/abcabc.cpp
    src/abcabc_verify.cpp
)
target_include_directories(pcg PUBLIC include)
target_link_libraries(pcg PUBLIC OpenSSL::Crypto)
target_compile_options(pcg PRIVATE -Wall -Wextra -Wpedantic)
add_executable(pcg_tests tests/test_main.cpp)
target_link_libraries(pcg_tests PRIVATE pcg)
enable_testing()
add_test(NAME pcg_tests COMMAND pcg_tests)
add_executable(pcg_golden_tests tests/golden_test.cpp)
target_link_libraries(pcg_golden_tests PRIVATE pcg)
add_test(NAME pcg_golden_tests COMMAND pcg_golden_tests ${CMAKE_SOURCE_DIR}/corpus/abcabc-conformance-v1.csv)

add_executable(pcg_sign_rotation_tests tests/sign_rotation_test.cpp)
target_link_libraries(pcg_sign_rotation_tests PRIVATE pcg)
add_test(NAME pcg_sign_rotation_tests COMMAND pcg_sign_rotation_tests ${CMAKE_SOURCE_DIR}/corpus/abcabc-sign-rotation-v1.csv)
