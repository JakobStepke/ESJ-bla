cmake_minimum_required(VERSION 3.20)
project(ASC_bla)

set (CMAKE_CXX_STANDARD 17)

include_directories(src)

find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE PYBIND11_DIR)
list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_DIR}")
find_package(pybind11 CONFIG REQUIRED)

find_package(LAPACK REQUIRED)
add_executable (test_vector tests/test_vector.cc)
target_link_libraries (test_vector PUBLIC ${LAPACK_LIBRARIES})
target_sources (test_vector PUBLIC src/vector.h src/expression.h)



find_package(LAPACK REQUIRED)
add_executable (test_lapack tests/test_lapack.cc)
target_link_libraries (test_lapack PUBLIC ${LAPACK_LIBRARIES})

find_package(LAPACK REQUIRED)
add_executable (all_tests tests/tests.cpp)
target_link_libraries (all_tests PUBLIC ${LAPACK_LIBRARIES})

add_executable (speed_test tests/speed_test.cc)
target_sources (speed_test PUBLIC src/vector.h src/expression.h)
target_link_libraries (speed_test PUBLIC ${LAPACK_LIBRARIES})

pybind11_add_module(bla src/bind_bla.cpp)

install (TARGETS bla DESTINATION ASCsoft)
install (FILES src/vector.h DESTINATION ASCsoft/include)
