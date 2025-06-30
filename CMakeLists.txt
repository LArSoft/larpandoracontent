cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

find_package(cetmodules 3.13.00 QUIET)

if(cetmodules_FOUND)
    set(LAR_PROJECT_NAME larpandoracontent)
else()
    set(LAR_PROJECT_NAME LArContent)
endif()

project(${LAR_PROJECT_NAME} VERSION 04.19.05 LANGUAGES CXX)

if(cetmodules_FOUND)
    message(STATUS "Configuring with cetmodules build system.")
    include(cmake/cetmodules_build.cmake)
else()
    message(STATUS "Configuring with standalone build system.")
    include(cmake/standalone_build.cmake)
endif()
