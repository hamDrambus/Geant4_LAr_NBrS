#!For finding Geant4 specify -DCMAKE_PREFIX_PATH=...geant4.10.04.p01-install/lib/Geant4-11.0.0
#!For finding Boost specify only root directory-DBOOST_ROOT=...Software/boost_1_67_0/
#!For finding Geant4 libraries during execution 'source .../geant4.10.04.p01-install/bin/geant4.sh'
#is required either at the start of each bash session or at system startup

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16)
project(Geant_simulation)

find_package(Boost 1.33 REQUIRED)
#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
#include(${Boost_USE_FILE}) does not work
#include(${Geant4_USE_FILE})
include_directories(${Boost_INCLUDE_DIR} ${Geant4_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/source/*.cc)
file (GLOB sources1 ${PROJECT_SOURCE_DIR}/source/*.cpp)
list (APPEND sources ${sources1})
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
file (GLOB headers1 ${PROJECT_SOURCE_DIR}/include/*.h)
list (APPEND headers ${headers1})

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(${PROJECT_NAME} ${sources} ${headers})
target_link_libraries(${PROJECT_NAME} PUBLIC ${Geant4_LIBRARIES} ${Boost_LIBRARIES})
message(STATUS "Geant4 libraries: ${Geant4_LIBRARIES}")

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build project. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(PROJECT_SCRIPTS
  vis.mac
  )

foreach(_script ${PROJECT_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${CMAKE_BUILD_TYPE}/${_script}
    COPYONLY
    )
endforeach()

add_definitions(${Geant4_DEFINITIONS} -DG4VIS_USE)
set(CMAKE_CXX_FLAGS ${Geant4_CXX_FLAGS})
set(CMAKE_CXX_STANDARD ${Geant4_CXXSTD})

#----------------------------------------------------------------------------
# Install the executable to Debug/Release directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_BINARY_DIR}/${CMAKE_BUILD_TYPE})
