cmake_minimum_required(VERSION 3.5)
project(basic_plus_ImageSeg LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

find_package(Qt6 COMPONENTS Core Gui Widgets PrintSupport REQUIRED)

set(SOURCE_FILES
    customgraphicsscene.cpp
    jetcolormap.cpp
    main.cpp
    ultrasound_cscan_seg.cpp
    ../basic_read/trimdialog.cpp
    ../basic_read/ultrasound_cscan_process.cpp
    ../basic_read/utils.cpp
    ../basic_read/qcustomplot.cpp
    ../basic_read/orthosliceviewer.cpp
)

set(HEADER_FILES
    customgraphicsscene.h
    jetcolormap.h
    ultrasound_cscan_seg.h
    ../basic_read/trimdialog.h
    ../basic_read/fftw3.h
    ../basic_read/ultrasound_cscan_process.h
    ../basic_read/utils.h
    ../basic_read/qcustomplot.h
    ../basic_read/orthosliceviewer.h
)


add_executable(${PROJECT_NAME} ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(${PROJECT_NAME} Qt6::Core Qt6::Gui Qt6::Widgets Qt6::PrintSupport)

# --------- FFTW3 linking if necessary
# CMakeLists.txt is in a directory at the same level as basic_read
set(FFTW_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../basic_read")
set(FFTW_LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../basic_read")
# Include FFTW3 header files
include_directories(${FFTW_INCLUDE_DIR})
# Directly specify the path to the FFTW3 library
set(FFTW_LIBRARY "${FFTW_LIB_DIR}/libfftw3-3.dll")
# Find FFTW3 Library
find_library(FFTW_LIBRARY NAMES fftw3 PATHS ${FFTW_LIB_DIR})
if(NOT FFTW_LIBRARY)
    message(FATAL_ERROR "FFTW3 library not found")
endif()
# Linking the FFTW3 library
target_link_libraries(${PROJECT_NAME} ${FFTW_LIBRARY})
# ----------

include(GNUInstallDirs)
install(TARGETS ${PROJECT_NAME}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
