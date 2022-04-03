Dependencies:
c++ BOOST
Geant4 (only v10.4 was used and tested)
CERN ROOT is required to analyze obtained data with root_scripts/* but it is optional.

Run cmake to create makefile:

mkdir Build
cd Build
cmake -DCMAKE_PREFIX_PATH=${HOME}/Software/Geant4/geant4.10.04.p01-install/lib/Geant4-11.0.0 ../
make

=================================================
Or run script to create Eclipse project:

bash GenEclipseProject.sh

Open in eclipse as File->Improt->General->Existing project in folder->browse to generated -build folder. Build via Project Explorer->Build Targets. Debug as C/C++ remote Application, with set binary location and disabled auto build. When any files or libraris are added to/revoved from the project, it must be regenerated with GenEclipseProject.sh.
