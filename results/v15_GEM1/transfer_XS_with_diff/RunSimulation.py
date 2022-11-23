#!/usr/bin/env python3
import sys
import os
import subprocess
from subprocess import PIPE
import shutil
import fileinput

current_path = os.path.dirname(os.path.realpath(__file__))
parent_path = os.path.dirname(current_path)
sys.path.append(parent_path)

import SolveFields

# These folders are relative to this script
ResultsFolder = "./"
SettingsTemplate = "./settings_template.xml"
Binary = "../../../../NBrS_THGEM_LAr_v0-build/RelWithDebInfo/Geant_simulation"
DataPath = "../../../data"
RecalculateField = False
RecalculateMesh = False
#V0s = [20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0]
Vt1s = [80, 90, 100, 110, 120, 130, 139, 160, 185, 200, 231, 250, 277, 300, 323, \
350, 369, 390, 416, 450, 462, 480, 508, 525, 554, 580, 600, 620, 646, 670, 693, \
620, 739, 760, 785, 800, 831, 850, 877, 900, 923, 950, 1000, 1050, 1100]
V0s = [20.0 for i in range(len(Vt1s))]

# returns binary absolute path, settings absolute path and log file absolute path (str1, str2, str3)
def prepare_settings(V0, Vth1):
    path = os.path.abspath(__file__)
    dir_path = os.path.dirname(path)
    abs_binary = os.path.normpath(os.path.join(dir_path, Binary))
    if not os.path.isfile(abs_binary):
        print ("ERROR: \""+abs_binary+"\" binary does not exist. Cannot run simulation, aborting.")
        return None
    abs_binary_path = os.path.dirname(abs_binary)
    abs_set_templ = os.path.normpath(os.path.join(dir_path, SettingsTemplate))
    if not os.path.isfile(abs_set_templ):
        print ("ERROR: \""+abs_set_templ+"\" settings template does not exist. Cannot set up simulation, aborting.")
        return None

    global RecalculateMesh
    field_map_files = SolveFields.solve_fields(V0, Vth1, recalculateField=RecalculateField, recalculateMesh=RecalculateMesh)
    if field_map_files is None:
        return None
    RecalculateMesh = False # Mesh needs to be recalulated only once

    suffix = str(round(Vth1)) + "V"
    abs_output_path = os.path.normpath(os.path.join(os.path.join(dir_path, ResultsFolder), suffix))
    os.makedirs(abs_output_path, exist_ok=True)
    abs_settings = os.path.normpath(os.path.join(os.path.join(dir_path, ResultsFolder), "settings_" + suffix + ".xml"))
    shutil.copy2(abs_set_templ, abs_settings)
    abs_data_path = os.path.normpath(os.path.join(dir_path, DataPath))
    rel_data_path = os.path.relpath(abs_data_path, abs_binary_path) # data path/folder is relative to folder from with binary is to be executed!
    rel_data_path = rel_data_path + "/" # Folders are stored with '/' in my c++ program
    rel_output_path = os.path.relpath(abs_output_path, abs_binary_path) # output path/folder is relative to folder from which binary is to be executed!
    rel_output_path = rel_output_path + "/"

    rel_mesh_path = os.path.relpath(field_map_files[0], abs_data_path) # The rest of folders in settings are relative to data folder
    rel_mesh_path = rel_mesh_path + "/"
    rel_field_map_file = os.path.relpath(field_map_files[1], abs_data_path)

    for line in fileinput.input([abs_settings], inplace=True):
        l = line.replace('DATA_PATH', rel_data_path) # binary is executed from its folder
        l = l.replace('OUTPUT_FOLDER', rel_output_path)
        l = l.replace('MESH_FOLDER', rel_mesh_path)
        print(l.replace('FIELD_MAP_FILE', rel_field_map_file), end='')
    abs_logfile = os.path.join(abs_output_path, "Log.txt")
    return (abs_binary, abs_settings, abs_logfile)


if __name__ == "__main__":
    if len(V0s) != len(Vt1s):
        print("ERROR: V0s and Vt1s have diffrent lengths (" + len(V0s) + " vs " + len(Vt1s) + "). Aborting RunSimulation.py.")
        sys.exit(1)
    for i in range(min(len(V0s), len(Vt1s))):
        print ("****************************************************************")
        print ("****************************************************************")
        print ("Starting simulating V0 = " + str(round(V0s[i], 1)) + ", Vthgem = " + str(round(Vt1s[i])) )
        print ("****************************************************************")
        print ("****************************************************************")
        files = prepare_settings(V0s[i], Vt1s[i])
        if files is None:
            continue
        abs_binary_path = os.path.dirname(files[0])
        binary = os.path.basename(files[0])
        binary = os.path.join("./", binary)
        """
        p1 = subprocess.Popen([binary, files[1]], stdout=PIPE, stderr=subprocess.STDOUT, cwd=abs_binary_path)
        p2 = subprocess.Popen(["tee", files[2]], stdin=p1.stdout)
        p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
        p2.communicate()
        if not p1.returncode is None:
            print ("****************************************************************")
            print ("****************************************************************")
            print("ERROR while running " + binary + ". Skipping V0= " + str(round(V0s[i], 1)) + ", Vthgem= " + str(round(Vt1s[i])) )
            print ("****************************************************************")
            print ("****************************************************************")
            continue
        """
        print ("****************************************************************")
        print ("****************************************************************")
        print ("Ended simulating V0 = " + str(round(V0s[i], 1)) + ", Vthgem = " + str(round(Vt1s[i])) )
        print ("****************************************************************")
        print ("****************************************************************")
    print("End of RunSimulation.py")
