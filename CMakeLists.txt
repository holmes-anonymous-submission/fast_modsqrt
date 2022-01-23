cmake_minimum_required(VERSION 3.17)
project(fast_modsqrt CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(GMP REQUIRED)
include_directories(${GMP_INCLUDE_DIR})
find_package(FLINT REQUIRED)
include_directories(${FLINT_INCLUDE_DIR}/flint)
find_package(MPFR REQUIRED)
include_directories(${MPFR_INCLUDE_DIR})

find_package(OpenMP)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")
else()
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -fopenmp")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

add_executable(fast_modsqrt main.cpp)
target_link_libraries(fast_modsqrt ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${FLINT_LIBRARIES} OpenMP::OpenMP_CXX)