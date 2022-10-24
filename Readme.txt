Dependencies:
1) c++ Boost. Headers are enough, no need for compiling it.
2) Geant4 (only v11.0 was used and tested). Do not forget to run '$source ../bin/geant4.sh' and to add Geant4's library path to LD_LIBRARY_PATH. Best to add 'source ../bin/geant4.sh' to startup script in /etc/profile.d/ 

Required for the project overall but not used in compilation of C++ code:
1) CERN ROOT is required to analyze obtained data with root_scripts/*.
2) FreeCAD, view3dscene or something else is recommended to view VRML2FILE (g4.wrl) visualization output. Otherwise refer to Geant4 visualization (configured from vis.mac)
3) Gmsh v3 is required to create THGEM cell mesh (.geo file->.msh)
4) Elmer is required to calculate electric field in cell using Gmsh output (.msh + .sif -> .results, .header, .nodes, .elements, .boundary).
5) gnuplot is recommended to quickly plot some simple files (such as electric fields). There are some funcitons which plot data by connecting to gnuplot with pipe. Not required for running the application.

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
=================================================
How to run simulation

Generally speaking simulation is configured from several places:
1) GlobalParamters (gPars::) is responsible for straightforward program setups such as filenames, some detector dimesnsions (which can't be fully arbitrary), number of events to simulate, and some debug/output options. UPD: Detector demensions and some global parameters related to them are now stored in DetectorSettings class.
2) UserInitialization is responsible for selection of detector geometry (there are the main one and a few for testing) and how and what particles are generated.
3) ArgonPropertiesTables is a hard-coded computation of several LAr paremters within certain boundaries. ArgonPropertiesTables::Initialize contains initialization of integration intervals and step sizes. These should be carefully changed if different input cross-sections are used (e.g. another gas) or paramters are needed for different electric fields or at different accuracy. Plotting of paramters also contains some hard-coded paramters (mainly axes ranges). Realistically, this class is the most challenging to configure correctly, but normally this should be done once and then forgotten about.
4) Electric field maps are input data which are calculated by 3-rd party Gmsh v3 and Elmer programs. Their input files are present in the project, but output must be generated separately (it is quite memory-heavy) before running simulation.
5) Each VDetectorContruction inheritor has SetSizeAndPosition(), Construct() and other virtual methods which contain construction of specific geometry. This is usual Geant4 business. Note that some GlobalParamters (DetectorSettings) are changed there so that primary particle generators work correcly for any any geometry without any adjustments in them or GlobalParamters.
6) VGeneratePrimaries inheritors are best configured by specifying Pattern from UserInitialization if possible. Otherwise, either add another pattern, create new class or change code inside existing ones. I.e. avoid cryptic hard-coding inside UserInitialization.
7) There is vis.mac which configures Geant4 visualization. Refer to geant4's docs.

1) and 2) are easiest to change without breaking anything and are supposed to be chagned most often.
3) is requried if physical properties must be changed.
4) is requried for field maps. Must be done for newly loaded project.
5) is for new geometries (new detector versions).
6) is primarily for testing. Or for some completely new simulation.

UPD: parameters are now set in setting.xml file
run as:
...build/RelWithDebInfo/Geant_simulation path/to/settings.xml | tee path/to/log.txt

UPD: There are now bash scripts streamlining simulation process.
One script calculates electric fields for several THGEM voltages
Another one then generates several settings files from template that will use those electric fields and run simulation program with the settings.

UPD: results/v11 and results/v10 now have RunSimulation.py which will automatically build mesh if necessary, calculate electric fields and then run geant4 simulation. 
