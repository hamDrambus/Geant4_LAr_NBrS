#!/bin/bash
#Script preparing Elmer's .sif file based on input voltage
#Input must be THGEM1 divider voltage
#Returns full path to prepared .sif file

ElmerFolder=Elmer_v00.02
#Detector constants:
THGEM1_Rtotal=$(echo "scale=4; 1.0+8.6+1.2" | bc)
THGEM1_Rtop=$(echo "scale=4; 1.0+8.6" | bc)
THGEM1_Rbot=$(echo "scale=4; 1.0" | bc)
V0_kV=-20.0
THGEM0_DIVIDER=$(echo "scale=4; 600.0/800.0" | bc)
#distances from THGEM1 electrode to THGEM0 and anode in real detector
EL_GAP=2.20
ANODE_THGEM1_dz=0.5
#distances from THGEM1 electrode to cathode and anode in Gmsh model
CATHODE_dZ=0.3
ANODE_dZ=0.3

if [ -z "$1" ];
then
	echo "Error: no parameter passed."
	exit 1
fi
V=$1
[ "$V" -eq "$V" ] 2>/dev/null
if [ $? -ne 0 ]; then
  echo "Error: paramter $var is not a number"
	exit 2
fi

V_top_THGEM1=$(echo "scale=4; ${V}*${THGEM1_Rtop}/${THGEM1_Rtotal}" | bc)
V_bot_THGEM1=$(echo "scale=4; ${V}*${THGEM1_Rbot}/${THGEM1_Rtotal}" | bc)
Edrift=$(echo "scale=4; (${V_bot_THGEM1} - ${V0_kV}*1000*${THGEM0_DIVIDER})/${EL_GAP}" | bc)
Einduction=$(echo "scale=4; (0 - ${V_top_THGEM1})/${ANODE_THGEM1_dz}" | bc)
Vcathode=$(echo "scale=4; ${V_bot_THGEM1} - ${Edrift}*${CATHODE_dZ}" | bc)
Vanode=$(echo "scale=4; ${V_top_THGEM1} + ${Einduction}*${ANODE_dZ}" | bc)

# Absolute path to this script, e.g. /home/user/bin/foo.sh
SCRIPT=$(readlink -f "$0")
# Absolute path this script is in, thus /home/user/bin
SCRIPTPATH=$(dirname "$SCRIPT")
OUTPATH=${SCRIPTPATH}/${ElmerFolder}
TEPMFILE=${SCRIPTPATH}/case_template.sif
VERSION=${V}V
OUTFILE=${OUTPATH}/case_${VERSION}.sif
cp $TEPMFILE $OUTFILE

#Find string and replace in file
sed -i "s/CASE_VERSION/$VERSION/" $OUTFILE
sed -i "s/CATHODE_POTENTIAL/$Vcathode/" $OUTFILE
sed -i "s/ANODE_POTENTIAL/$Vanode/" $OUTFILE
sed -i "s/TOP_THGEM_POTENTIAL/$V_top_THGEM1/" $OUTFILE
sed -i "s/BOT_THGEM_POTENTIAL/$V_bot_THGEM1/" $OUTFILE

echo $OUTFILE
