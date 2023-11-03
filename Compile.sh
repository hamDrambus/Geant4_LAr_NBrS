#!/bin/sh

# Necessary input for CMake
GEANT4PATH=${HOME}/Software/Geant4/geant4-v11.0.0-install/lib/Geant4-11.0.0/
BOOSTPATH=${HOME}/Software/boost_1_67_0/

# Absolute path to this script, e.g. /home/user/bin/foo.sh
SCRIPT=$(readlink -f "$0")
# Absolute path this script is in, thus /home/user/bin
SCRIPTPATH=$(dirname "$SCRIPT")
PROJECTNAME=$(basename "$SCRIPTPATH")

# Clear build directory and cd to it
rm -rf ${SCRIPTPATH}-build
mkdir ${SCRIPTPATH}-build
cd ${SCRIPTPATH}-build

#BUILDTYPE=Debug
BUILDTYPE=RelWithDebInfo
#BUILDTYPE=Release

# set -x displays cmake command. Brackets create subshell so that there is no need ot call set +x
(set -x; cmake -DCMAKE_BUILD_TYPE=${BUILDTYPE} -DCMAKE_PREFIX_PATH=${GEANT4PATH} -DBOOST_ROOT=${BOOSTPATH} ../${PROJECTNAME})
(set -x; make -j6 install)


