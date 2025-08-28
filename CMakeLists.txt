cmake_minimum_required(VERSION 3.5)

project(DDMethod LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(GNUInstallDirs)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    add_definitions(-DUSE_OPENMP)
endif()

set(SOURCES
    DDNumericalMethod.cpp
)

set(HEADERS
    DDNumericalMethod.h
    GausLegendreQuadrature.h
    LegendrePolynomial.h
    BuildMatricesWIthTxt.h
    VariablesUsed.h
    get_GQ.h
)

add_library(DDMethod SHARED ${SOURCES} ${HEADERS})

target_include_directories(DDMethod
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(DDMethod PUBLIC OpenMP::OpenMP_CXX atomic)
else()
    target_link_libraries(DDMethod PUBLIC atomic)
endif()

install(TARGETS DDMethod
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Executável de teste (main.cpp)
add_executable(DDMethodTest main.cpp)
target_link_libraries(DDMethodTest PRIVATE DDMethod)
