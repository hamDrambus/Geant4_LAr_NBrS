Dependencies:
1) c++ Boost. Headers are enough, no need for compiling it.
2) Geant4 (only v10.4 was used and tested). Do not forget to run '$source ../bin/geant4.sh' and to add Geant4's library path to LD_LIBRARY_PATH. Best to add 'source ../bin/geant4.sh' to startup script in /etc/profile.d/ and 

Required for the project overall but not used in compilation of C++ code:
1) CERN ROOT is required to analyze obtained data with root_scripts/*.
2) FreeCAD or something else is required to view VRML2FILE (g4.wrl) visualization output.
3) Gmsh v3 is required to create THGEM cell mesh (.geo file->.msh)
4) Elmer is required to calculate electric field in cell using Gmsh output (.msh + .sif -> .results, .header, .nodes, .elements, .boundary).
5) gnuplot is recommended to quickly plot some simple files (such as electric fields).

=================================================
To compile Geant4 code run cmake to create makefile:

mkdir Build
cd Build
cmake -DCMAKE_PREFIX_PATH=${HOME}/Software/Geant4/geant4.10.04.p01-install/lib/Geant4-11.0.0 ../
make

Or run script to create Eclipse project:

bash GenEclipseProject.sh

Open in eclipse as File->Improt->General->Existing project in folder->browse to generated -build folder. Build via Project Explorer->Build Targets. Debug as C/C++ remote Application, with set binary location and disabled auto build. When any files or libraris are added to/revoved from the project, it must be regenerated with GenEclipseProject.sh.

=================================================
Note that Eclipse console is somewhat buggy. Because of that, progress bars behave incorrectly. To remedy this there are tow options:
1) debug application as a remote, i.e. launch it in native terminal and attach debugger to it from eclipse after. This requires pause at the start of the program so it not optimal.
2) Another option is to configure gbd (debugger) and and launcher to use native terminal instead. On ubuntu this is done using tty devices (emulated as text files). In Debug Configuration->Common set Input and Output files to open terminal's device file (e.g. /dev/pts/1), which can be obtained from open terminal by running 'tty'. In this case this terminal will be connected to and intercepted by Eclipse.
