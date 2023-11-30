cmake_minimum_required(VERSION 3.5)

project(basic_read LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(Qt6 COMPONENTS Core Gui Widgets PrintSupport REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# FFTW3 Header and Library Path
set(FFTW_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(FFTW_LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
# Directly specify the path to the FFTW3 library
set(FFTW_LIBRARY "${FFTW_LIB_DIR}/libfftw3-3.dll")
# Find FFTW3 Library
find_library(FFTW_LIBRARY NAMES fftw3 PATHS ${FFTW_LIB_DIR})
if(NOT FFTW_LIBRARY)
    message(FATAL_ERROR "FFTW3 library not found")
endif()

add_executable(basic_read
    main.cpp
    npy.hpp
    qcustomplot.cpp
    qcustomplot.h
    trimdialog.cpp
    trimdialog.h
    ultrasound_cscan_process.cpp
    ultrasound_cscan_process.h
    utils.cpp
    utils.h)

# Include FFTW3 Headers
target_include_directories(basic_read PRIVATE ${FFTW_INCLUDE_DIR})

# Linking Libraries
target_link_libraries(basic_read ${FFTW_LIBRARY})
target_link_libraries(basic_read Qt6::Core Qt6::Gui Qt6::Widgets Qt6::Widgets Qt6::PrintSupport)

include(GNUInstallDirs)
install(TARGETS basic_read
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
