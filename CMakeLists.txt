cmake_minimum_required(VERSION 3.16)
if(POLICY CMP0144)
  cmake_policy(SET CMP0144 NEW)
endif()
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
project(vina_aco CXX)

# Vina-ACO: AutoDock Vina with an Ant Colony Optimization global search (--search aco).
# One-command build: ./build.sh   (outputs bin/vina_aco)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
option(VINA_ACO_NATIVE "Compile with -march=native" ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
if(VINA_ACO_NATIVE)
  add_compile_options(-march=native)
endif()

set(Boost_NO_WARN_NEW_VERSIONS ON)
find_package(Boost 1.60 REQUIRED COMPONENTS program_options filesystem system thread serialization)
find_package(Threads REQUIRED)

set(VINA_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src)
file(GLOB VINA_LIB_SOURCES ${VINA_SRC}/lib/*.cpp)

add_library(vina_aco_lib STATIC ${VINA_LIB_SOURCES})
target_include_directories(vina_aco_lib PUBLIC ${VINA_SRC}/lib)
target_include_directories(vina_aco_lib SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
target_compile_definitions(vina_aco_lib PUBLIC VERSION="v1.2.7-aco")
target_link_libraries(vina_aco_lib PUBLIC Boost::program_options Boost::filesystem Boost::system
                      Boost::thread Boost::serialization Threads::Threads)

add_executable(vina_aco ${VINA_SRC}/main/main.cpp)
target_link_libraries(vina_aco PRIVATE vina_aco_lib)
set(VINA_ACO_BIN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/bin" CACHE PATH "Where the vina_aco executable is written")
set_target_properties(vina_aco PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${VINA_ACO_BIN_DIR}
  BUILD_RPATH "${Boost_LIBRARY_DIRS}"
  INSTALL_RPATH "${Boost_LIBRARY_DIRS}")
install(TARGETS vina_aco RUNTIME DESTINATION bin)

# Math invariants of the ACO internals (tests/test_aco_math.cpp includes aco.cpp on purpose;
# assertions are kept on even in Release builds).
option(VINA_ACO_TESTS "Build the ACO math test" ON)
if(VINA_ACO_TESTS AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_aco_math.cpp)
  enable_testing()
  add_executable(test_aco_math tests/test_aco_math.cpp)
  target_compile_options(test_aco_math PRIVATE -UNDEBUG)
  target_link_libraries(test_aco_math PRIVATE vina_aco_lib)
  set_target_properties(test_aco_math PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                        BUILD_RPATH "${Boost_LIBRARY_DIRS}")
  add_test(NAME aco_math COMMAND test_aco_math)
endif()
