#!/bin/bash

#folders are relative to this script
ResultsFolder=results/v8.1_no_diff
Binary=../NBrS_THGEM_LAr_v0-build/RelWithDebInfo/Geant_simulation
declare -a Voltages=(6180 5993 5728 5297 4856 4413 3972 3531 3090 2648 2206 1765)

	# Absolute path to this script, e.g. /home/user/bin/foo.sh
SCRIPT=$(readlink -f "$0")
# Absolute path this script is in, thus /home/user/bin
SCRIPTPATH=$(dirname "$SCRIPT")
OUTPATH=${SCRIPTPATH}/${ResultsFolder}

function prepare_settings () {
# $1 is voltage (e.g. string 6180)
local TEPMFILE=${OUTPATH}/settings_template.xml
local OUTFILE=${OUTPATH}/settings_${1}V.xml
cp $TEPMFILE $OUTFILE

#Find string and replace in file
sed -i "s/VOLTAGE/${1}/" $OUTFILE
mkdir -p ${OUTPATH}/${1}V/
echo ${OUTFILE}
}

if [ -z ${ResultsFolder+x} ]; then
  echo "ResultsFolder is unset"
  exit 1
else
  echo "ResultsFolder is set to '$ResultsFolder'"
fi

if [ -z ${Binary+x} ]; then
  echo "Binary is unset"
  exit 1
else
  echo "Binary is set to '$Binary'"
fi

cd $(dirname "$Binary")
Bin=$(basename "$Binary")

for V in ${Voltages[@]}; do
  Settings=$(prepare_settings $V)
  if [ -z ${Settings+x} ]; then
	echo "Error: Could not prepare settings for ${V} volts"
  else
    Log=${OUTPATH}/${V}V/Log.txt
    LogPath=$(dirname "${Log}")
    mkdir -p ${LogPath}
    echo "Executing \$(./${Bin} ${Settings} | tee ${Log})"
    (./${Bin} ${Settings} | tee ${Log})
  fi
done

