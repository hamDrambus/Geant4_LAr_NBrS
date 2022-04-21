#!/bin/bash

declare -a Voltages=(6180 5993 5728 5297 4856 4413 3972 3531 3090 2648 2206 1765)
declare -a SifNames=()
for V in ${Voltages[@]}; do
  sifFile=$(./PrepareSif.sh $V)
  if [ -z ${ElmerFolder+x} ]; then
    ElmerFolder=$(dirname "$sifFile")
  fi
  SifNames+=( $(basename "$sifFile") )
done

if [ -z ${ElmerFolder+x} ]; then
  echo "ElmerFolder is unset"
  exit 1
else
  echo "ElmerFolder is set to '$ElmerFolder'"
fi

cd ${ElmerFolder}
for sif in ${SifNames[@]}; do
  echo "Executing \$(ElmerSolver $sif)"
  (ElmerSolver $sif | tee Log.txt)
done
