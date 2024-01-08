cmake_minimum_required(VERSION 3.5)

project(DDMethod LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Lista de arquivos de cabeçalho
set(HEADERS
    DataMatrices.h
    DDNumericalMethold.h
    GausLegendreQuadrature.h
    LegendrePolynomial.h
    VariablesUsed.h
)

# Lista de arquivos de origem
set(SOURCES
    main.cpp
)

# Criar a biblioteca compartilhada (.so)
add_library(DDMethod SHARED ${HEADERS} ${SOURCES})

# Adicionar os diretórios de inclusão
target_include_directories(DDMethod
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

# Instalar a biblioteca compartilhada no diretório especificado
install(TARGETS DDMethod
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
