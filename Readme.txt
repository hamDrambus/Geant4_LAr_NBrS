Dependencies:
1) c++ Boost. Headers are enough, no need for compiling it.
2) Geant4 (only v11.0 was used and tested). Do not forget to run '$source ../bin/geant4.sh' and to add Geant4's library path to LD_LIBRARY_PATH. Best to add 'source ../bin/geant4.sh' to startup script in /etc/profile.d/ 

Required for the project overall but not used in compilation of C++ code:
1) CERN ROOT is strongly recommended to analyze obtained data with root_scripts/*. Alternatively, it can be done by hand using text output of Geant4_simulation (do note QE and PDE issue discussed below).
2) FreeCAD, view3dscene or something else is recommended to view VRML2FILE (g4.wrl) visualization output. Otherwise refer to Geant4 visualization (configured from vis.mac)
3) Gmsh v4 is required to create THGEM cell mesh (.geo file->.msh). Gmsh v3 does not have mesh options setups from .geo file (could not find its documentation) so it won't work properly.
4) Elmer is required to calculate electric fields in cell using Gmsh output (.msh + .sif -> .results, .header, .nodes, .elements, .boundary).
5) Python3 is strongly recommended to run simulation scripts (RunSimulation.py) which streamline the workflow. In particular, they invoke Gmsh and Elmer if correct script set-up and input files are used. After electric fields are calculated, settings.xml is automatically generated from template settings using the field maps files. Finally, c++ geant4 simulation is run using generated field maps and settings. The scripts are very readable and thus new simulation cases can be added as needed.
6) gnuplot is recommended to quickly plot some simple files (such as electric fields). There are some functions in the c++ program which plot data by connecting to gnuplot with pipe. Not required for running the application.

=================================================
To compile Geant4 code run cmake to create makefile:

mkdir build
cd build
cmake -DCMAKE_PREFIX_PATH=${HOME}/Software/Geant4/geant4.10.04.p01-install/lib/Geant4-11.0.0 ../
make
make install

Or run script to create Eclipse project:

bash GenEclipseProject.sh

Open in eclipse as [File->Import->General->Existing project in folder]->browse to generated -build folder. Build via [Project Explorer->Build Targets]. Debug as 'C/C++ Attach to Application' or 'C/C++ Application', with set binary location and disabled auto build (if necessary). When any files or libraries are added to/removed from the project, it must be regenerated with GenEclipseProject.sh.

=================================================
Note that Eclipse console is somewhat buggy. Because of that, progress bars behave incorrectly. To remedy this there are two options:
1) debug application as a remote, i.e. launch it in native terminal and attach debugger to it from eclipse after. This requires pause at the start of the program so it not optimal.
2) Another option is to configure gbd (debugger) and launcher to use native terminal instead. On ubuntu this is done using tty devices (emulated as text files). In [Debug Configuration->Common] set Input and Output files to open terminal's device file (e.g. /dev/pts/1), which can be obtained from open terminal by running 'tty'. In this case this terminal will be connected to and intercepted by Eclipse.

=================================================
How to run simulation from scratch (github):
	1) Install required and recommended dependencies.
	2) Load the code and compile it using CMake or Eclipse.
	4) Check paths to binaries (Geant4_simulation, Elmer, Gmsh) and folders in python scripts SolveFields.py and RunSimulation.py. These are set up at the head of the scripts.
	5) python3 project_folder/results/v*/*/RunSimulation.py
	6) Process results using project_folder/root_scripts/*. E.g:
		cd project_folder/root_scripts
		root -l
		.L init.cpp
		.L plot_Npe_spectrum.cpp
		.x print_Npe_vs_V.cpp
	Alternatively, log files generated during execution of RunSimulation.py can be used directly by hand.

=================================================
Setting up the simulation

Generally speaking simulation is configured from several places:
1) GlobalParamters (gPars::) is responsible for straightforward program setups such as filenames, some detector dimesnsions (which can't be fully arbitrary), number of events to simulate, and some debug/output options. UPD: Detector demensions and some global parameters related to them are now stored in DetectorSettings class.
2) UserInitialization is responsible for selection of detector geometry (there are the main one and a few for testing) and how and what particles are generated.
3) ArgonPropertiesTables is a hard-coded computation of several LAr paremters within certain boundaries. ArgonPropertiesTables::Initialize contains initialization of integration intervals and step sizes. These should be carefully changed if different input cross-sections are used (e.g. another gas) or paramters are needed for different electric fields or at different accuracy. Plotting of paramters also contains some hard-coded paramters (mainly axes ranges). Realistically, this class is the most challenging to configure correctly, but normally this should be done once and then forgotten about.
4) Electric field maps are input data which are calculated by 3-rd party Gmsh and Elmer programs. Their input files are present in the project, but output must be generated separately (it is quite memory-heavy) before running simulation.
5) Each VDetectorContruction inheritor has SetSizeAndPosition(), Construct() and other virtual methods which contain construction of specific geometry. This is usual Geant4 business. Note that some GlobalParamters (DetectorSettings) are changed there so that primary particle generators work correctly for any geometry without any additional adjustments in them or GlobalParamters. Do also note that PMTs' quantum efficiency (QE) and SiPMs' photon detection efficiency (PDE) are by default set to 1 in geant4. It is done to increase number of photons and decrease statistical errors. QE and PDE are taken into account (averaged over detected/emitted spectrum) in root_scripts/*. 
6) VGeneratePrimaries inheritors are best configured by specifying Pattern from UserInitialization if possible. Otherwise, either add another pattern, create new class or change code inside existing ones. I.e. avoid cryptic hard-coding inside UserInitialization.
7) There is vis.mac which configures Geant4 visualization. Refer to geant4's docs.

1) and 2) are easiest to change without breaking anything and are supposed to be changed most often.
3) is required if physical properties must be changed.
4) is required for field maps. Must be adjusted/set for newly loaded project.
5) is for new geometries (new detector versions).
6) is primarily for testing. Or for some completely new simulation.

UPD: parameters are now set in settings.xml file.
Thus program must be run as:
...build/RelWithDebInfo/Geant_simulation path/to/settings.xml | tee path/to/log.txt
If no parameter is passed, then "settings.xml" is expected to be present in the current directory

UPD: There are now python scripts streamlining simulation process.
One script calculates electric fields for several THGEM voltages
Another one then generates several settings files from template that will use those electric fields and run simulation program with the settings.

UPD: results/v11 and results/v10 now have RunSimulation.py which will automatically build mesh if necessary, calculate electric fields and then run geant4 simulation. 

=================================================
Systematic errors legend:
Theory:
T1_R	Reflectance values: FR4, Cu, Steel, Wires
T2_D	Diffused reflections: Cu, FR4, LAr-gas, Steel wires, Acrilyc
T3_PDE	SiPM PDE uncertainty
T4_XY	Source position uncertainty
T5_XY	Source XY profile uncertainty
T6_Di	Electron diffusion uncertainty
T7_MSH	Electric field precision (mesh precision)
T8_XS	Cross-sections' uncertainties (e.g. discrepancy in drift velocity)
T9	Is stationary theory applicable in non-uniform field?

Experiment:
E1_Cal	SiPM calibration uncertainty
E2_Ev	Event selection
E3_St	Statistical error
E4_Qcal	Absolute values uncertainty (calibration for X-ray or theory for alphas)
E5_V0	Drift field uncertainty
E6_Q_S1	Charge instability by S1
E7_Qinst	Charge instability per events (first 1k events vs all events) (same as above?)
E8_O2	O2 impurity
E9*	THGEM/GEM charge-up effects
E10_V1*	THGEM voltage uncertainty

