#!/bin/bash

#folders are relative to this script
ResultsFolder=./
Binary=../../../../NBrS_THGEM_LAr_v0-build/RelWithDebInfo/Geant_simulation
#declare -a V0s=(22.0 21.5 21.0 20.5 20.0 19.5 19.0 18.5 18.0 17.5 17.0 16.5 16.0 15.5 15.0 14.5 14.0 13.0 12.0 11.0)
#declare -a Vt1s=(0 0 0 0 0 0 0 0 338 563 900 1238 1238 1238 1238 1238 1688 2250 2250 2025)

declare -a V0s=(11.0)
declare -a Vt1s=(2025)

# Absolute path to this script, e.g. /home/user/bin/foo.sh
SCRIPT=$(readlink -f "$0")
# Absolute path this script is in, thus /home/user/bin
SCRIPTPATH=$(dirname "$SCRIPT")
OUTPATH=${SCRIPTPATH}/${ResultsFolder}

function prepare_settings () {
  # $1 is V0 voltage (e.g. string 22.0) and $2 is Vthgem1 voltage (on the divider)
  local TEPMFILE=${OUTPATH}/settings_template.xml
  local OUTFILE=${OUTPATH}/settings_${1}V.xml
  cp $TEPMFILE $OUTFILE

  #Find string and replace in file
  sed -i "s/VOLTAGE/${1}/" $OUTFILE
  mkdir -p ${OUTPATH}/${1}V/
  echo ${OUTFILE}
}

if [ -z ${ResultsFolder+x} ]; then
  echo "ResultsFolder is unset" > /dev/stderr
  exit 1
else
  echo "ResultsFolder is set to '$ResultsFolder'" > /dev/stderr
fi

if [ -z ${Binary+x} ]; then
  echo "Binary is unset" > /dev/stderr
  exit 1
else
  echo "Binary is set to '$Binary'" > /dev/stderr
fi

# Absolute path to this script, e.g. /home/user/bin/foo.sh
SCRIPT=$(readlink -f "$0")
# Absolute path this script is in, thus /home/user/bin
SCRIPTPATH=$(dirname "$SCRIPT")
cd $(dirname "$SCRIPTPATH/$Binary")
Bin=$(basename "$Binary")

for i in "${!V0s[@]}"; do
  Settings=$(prepare_settings ${V0s[i]} ${Vt1s[i]})
  echo "Executing \$(bash ${SCRIPTPATH}/../SolveFields.sh ${V0s[i]} ${Vt1s[i]})" > /dev/stderr
  (${SCRIPTPATH}/../SolveFields.sh ${V0s[i]} ${Vt1s[i]})
  if [ -z ${Settings+x} ]; then
	  echo "Error: Could not prepare settings for ${V0s[i]} kV" > /dev/stderr
  else
    Log=${OUTPATH}/${V0s[i]}V/Log.txt
    LogPath=$(dirname "${Log}")
    mkdir -p ${LogPath}
    echo "Executing \$(./${Bin} ${Settings} | tee ${Log})" > /dev/stderr
    (./${Bin} ${Settings} | tee ${Log})
  fi
done

