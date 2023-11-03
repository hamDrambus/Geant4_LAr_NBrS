#!For finding Geant4 specify -DCMAKE_PREFIX_PATH=...geant4.11.0.0-install/lib/Geant4-11.0.0
#!For finding Boost specify only root directory-DBOOST_ROOT=...Software/boost_1_67_0/
#!For finding Geant4 libraries during execution 'source .../geant4.11.0.0-install/bin/geant4.sh'
#is required either at the start of each bash session or at system startup

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16)
project(Geant_simulation)
cmake_policy(SET CMP0144 NEW)

find_package(Boost 1.33 REQUIRED)
#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 10.4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 10.4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
#include(${Boost_USE_FILE}) does not work
#include(${Geant4_USE_FILE})
include_directories(${Boost_INCLUDE_DIR} ${Geant4_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/externals/indicators/include ${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
file(GLOB_RECURSE sources ${PROJECT_SOURCE_DIR}/source/*.cc)
file (GLOB_RECURSE sources1 ${PROJECT_SOURCE_DIR}/source/*.cpp)
list (APPEND sources ${sources1})
file(GLOB_RECURSE headers ${PROJECT_SOURCE_DIR}/include/*.hh)
file(GLOB_RECURSE headers1 ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB_RECURSE headers2 ${PROJECT_SOURCE_DIR}/include/*.hpp)
list (APPEND headers ${headers1})
list (APPEND headers ${headers2})

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

#----------------------------------------------------------------------------
# Set -rpath for binary. Allows it to find libraries in non-stardard locations
# (e.g. not in /usr/) without using LD_LIBRARY_PATH
#
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

#----------------------------------------------------------------------------
# Remove some of the warnings for gcc as
# boost has quite a lot of these particular ones
#
if (MSVC)
    # warning level 4
    add_compile_options(/W4)
else()
	# additional warnings
	add_compile_options(-Wno-pedantic -Wno-shadow -Wno-reorder -Wno-unused-parameter)
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
  .gdbinit
  )
list(TRANSFORM PROJECT_SCRIPTS PREPEND ${PROJECT_SOURCE_DIR}/)

add_definitions(${Geant4_DEFINITIONS} -DG4VIS_USE)
set(CMAKE_CXX_FLAGS ${Geant4_CXX_FLAGS})
set(CMAKE_CXX_STANDARD ${Geant4_CXXSTD})

#----------------------------------------------------------------------------
# Install the executable to Debug/Release directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE})
install(FILES ${PROJECT_SCRIPTS} DESTINATION ${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE})
