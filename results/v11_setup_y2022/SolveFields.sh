#!/bin/bash
# Input is V0 voltage [kV] ($1) and Vthgem1 voltage [V] ($2)

BinaryGmsh=~/Software/Gmsh_v3/gmsh-3.0.6-Linux64/bin/gmsh
BinaryElmer=ElmerGrid
ElmerFolder=Elmer_v00.01
Geo=v00.01_THGEM0.geo

# Absolute path to this script, e.g. /home/user/bin/foo.sh
SCRIPT=$(readlink -f "$0")
# Absolute path this script is in, thus /home/user/bin
SCRIPTPATH=$(dirname "$SCRIPT")
Geo=${SCRIPTPATH}/${Geo}

#Detector constants:
THGEM1_Rtotal=$(echo "scale=4; 1.0+8.6+1.2" | bc)
THGEM1_Rthgem=$(echo "scale=4; 1.0+8.6" | bc)

R1=60
R2=80
R3=260
R4=240
Rtot=$(echo "scale=4; ${R1}+3*${R2}+${R3}+${R4}" | bc)

#distances from THGEM0 electrode to THGEM1 and cathode in real detector
EL_GAP_FULL=2.20
EL_GAP=0.60
DRIFT_L=4.8
#LAr dielectric constant
LAR_EPS=1.54
#distances from THGEM1 electrode to cathode and anode in Gmsh model
CATHODE_dZ=0.5
ANODE_dZ=0.5

function is_number () {
  if [ -z "$1" ];
  then
	  return 1
  fi
  if ! [[ "$1" =~ ^[+-]?[0-9]+([.][0-9]+)?$ ]]; then
    echo "Error: parameter $1 is not a number" > /dev/stderr
	  return 1
  fi
  return 0
}

#input is only mesh folder to check ($1)
function do_need_meshing () {
  if [ -f ${1}/mesh.boundary ]; then
    return 1
  fi
  if [ -f ${1}/mesh.elements ]; then
    return 1
  fi  
  if [ -f ${1}/mesh.header ]; then
    return 1
  fi
  if [ -f ${1}/mesh.nodes ]; then
    return 1
  fi
  return 0
}

#returns mesh folder
function create_mesh() {
  GeoFile=$(basename "$Geo")
  GeoFolder=${SCRIPTPATH}/${GeoFile%.*}
  if do_need_meshing ${GeoFolder} ; then
    echo "Executing \$($BinaryGmsh $Geo -3 -order 2 -format msh)" > ${SCRIPTPATH}/Log_gmsh.txt
    ($BinaryGmsh $Geo -3 -order 2 -format msh > ${SCRIPTPATH}/Log_gmsh.txt)
    echo "Executing \$($BinaryElmer 14 2 ${GeoFolder}.msh -autoclean)" > ${SCRIPTPATH}/Log_gmsh.txt
    ($BinaryElmer 14 2 ${GeoFolder}.msh -autoclean > ${SCRIPTPATH}/Log_gmsh.txt)
    rm ${GeoFolder}.msh
  fi
  if do_need_meshing ${GeoFolder} ; then
    echo "Error: Failed to create mesh files." > ${SCRIPTPATH}/Log_gmsh.txt
	  exit 2
  fi  
  cp ${SCRIPTPATH}/diels.dat ${GeoFolder}/diels.dat
  echo ${GeoFolder}
}

if ! is_number "$1" ; then
	echo "Error: V0 (parameter #1) is invalid." > /dev/stderr
	exit 1
fi
if ! is_number "$2" ; then
	echo "Error: Vthgem1 (parameter #2) is invalid." > /dev/stderr
	exit 1
fi

V0=$1
Vt1=$2

MeshFolder=$(create_mesh)
echo "MeshFolder=${MeshFolder}"

#Calculating potentials for sif file
V_top_THGEM0=$(echo "scale=4; 1000*${V0}*${R3}/${Rtot}" | bc)
V_bot_THGEM0=0
Edrift=$(echo "scale=4; (${V0}*1000*3*${R2}/${Rtot})/${DRIFT_L}" | bc)
#Vgap in the real detector [V]
Vgap=$(echo "scale=4; (${V0}*1000*${R4}/${Rtot} + ${Vt1}*${THGEM1_Rthgem}/${THGEM1_Rtotal})" | bc)
Einduction=$(echo "scale=4; ${Vgap}/(${LAR_EPS}*${EL_GAP} + ${EL_GAP_FULL} - ${EL_GAP})" | bc)
Vcathode=$(echo "scale=4; ${V_bot_THGEM0} - ${Edrift}*${CATHODE_dZ}" | bc)
Vanode=$(echo "scale=4; ${V_top_THGEM0} + ${Einduction}*${ANODE_dZ}" | bc)

OUTPATH=${SCRIPTPATH}/${ElmerFolder}
TEPMFILE=${SCRIPTPATH}/case_template.sif
VERSION=${V0}v
SIFFILE=${OUTPATH}/case_${VERSION}.sif

mkdir -p $(dirname "$SIFFILE")
cp $TEPMFILE $SIFFILE

#Find string and replace in file
sed -i "s+MESH_FOLDER+$MeshFolder+" $SIFFILE
sed -i "s/CASE_VERSION/$VERSION/" $SIFFILE
sed -i "s/CATHODE_POTENTIAL/$Vcathode/" $SIFFILE
sed -i "s/ANODE_POTENTIAL/$Vanode/" $SIFFILE
sed -i "s/TOP_THGEM_POTENTIAL/$V_top_THGEM0/" $SIFFILE
sed -i "s/BOT_THGEM_POTENTIAL/$V_bot_THGEM0/" $SIFFILE

# Now is full path
ElmerFolder=$(dirname "$SIFFILE")

cd ${ElmerFolder}
echo "Executing \$(ElmerSolver $SIFFILE)" > /dev/stderr
(ElmerSolver $SIFFILE > ${ElmerFolder}/Log_elmer_${VERSION}.txt)

