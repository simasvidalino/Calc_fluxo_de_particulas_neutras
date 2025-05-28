cmake_minimum_required(VERSION 3.5)

project(DDMethod LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Incluir o módulo GNUInstallDirs
include(GNUInstallDirs)

# Lista de fontes para a biblioteca
set(SOURCES
    DDNumericalMethod.cpp
)

# Lista de cabeçalhos
set(HEADERS
    DDNumericalMethod.h
    GausLegendreQuadrature.h
    LegendrePolynomial.h
    BuildMatricesWIthTxt.h
    VariablesUsed.h
    get_GQ.h

)

# Crie a biblioteca compartilhada (.so)
add_library(DDMethod SHARED ${SOURCES} ${HEADERS})

# Inclua os diretórios de cabeçalho
target_include_directories(DDMethod
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Instale a biblioteca e os cabeçalhos
install(TARGETS DDMethod
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Crie o executável de teste
add_executable(DDMethodTest main.cpp)

# Link o executável com a biblioteca DDMethod
target_link_libraries(DDMethodTest PRIVATE DDMethod)
