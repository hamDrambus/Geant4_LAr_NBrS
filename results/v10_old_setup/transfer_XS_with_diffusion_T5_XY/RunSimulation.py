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
Vt1s = [200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, \
1600, 1900, 2100, 2400, 2800, 3150, 3300, 3700, 4200, 4600, 5050, 5500, 6400, 6600, \
6800, 7000, 7400, 7800, 8200, 8600, 9200]
Vt1s = Vt1s + [1765, 2206, 2648, 3090, 3531, 3972, 4413, 4856, 5297, 5738, 5993, 6180] # Voltages for 220203 experiment
Vt1s = Vt1s + [1762, 2644, 3525, 3966, 4847, 5288, 5728, 6169] # Voltages for 220804 experiment
#Vt1s = [6169]
Vt1s = sorted(list(set(Vt1s)))
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

        p1 = subprocess.Popen([binary, files[1]], stdout=PIPE, stderr=subprocess.STDOUT, cwd=abs_binary_path)
        p2 = subprocess.Popen(["tee", files[2]], stdin=p1.stdout)
        p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
        p2.communicate()
        os.remove(files[1])
        if not p1.returncode is None:
            print ("****************************************************************")
            print ("****************************************************************")
            print("ERROR while running " + binary + ". Skipping V0= " + str(round(V0s[i], 1)) + ", Vthgem= " + str(round(Vt1s[i])) )
            print ("****************************************************************")
            print ("****************************************************************")
            continue
        print ("****************************************************************")
        print ("****************************************************************")
        print ("Ended simulating V0 = " + str(round(V0s[i], 1)) + ", Vthgem = " + str(round(Vt1s[i])) )
        print ("****************************************************************")
        print ("****************************************************************")
    print("End of RunSimulation.py")
