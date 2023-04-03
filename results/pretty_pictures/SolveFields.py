#!/usr/bin/env python3
import os
import sys
import subprocess
import shutil
import fileinput

BinaryGmsh = "~/Software/Gmsh_v4.11.0/bin/gmsh"
BinaryElmer = "ElmerGrid"
BinarySolver = "ElmerSolver"
# Folder and file are relative to this script
ElmerFolder = "Elmer_v00.01"
GeoFile = "v00.01_THGEM1.geo"

def do_need_meshing (path):
    file = os.path.join(path, "mesh.boundary")
    if not os.path.isfile(file):
        return True
    file = os.path.join(path, "mesh.elements")
    if not os.path.isfile(file):
        return True
    file = os.path.join(path, "mesh.header")
    if not os.path.isfile(file):
        return True
    file = os.path.join(path, "mesh.nodes")
    if not os.path.isfile(file):
        return True
    return False

# return folder with geometry
def create_mesh(abs_geo_file, recalculate=False):
    abs_binary_gmsh = os.path.expanduser(BinaryGmsh)
    geo_filename = os.path.splitext(os.path.basename(abs_geo_file))[0]
    directory = os.path.dirname(abs_geo_file)
    abs_mesh_folder = os.path.join(directory, geo_filename)
    if do_need_meshing(abs_mesh_folder) or recalculate:
        if not os.path.isfile(abs_geo_file):
            print ("ERROR: \""+abs_geo_file+"\" file does not exist")
            return None
        print("Executing \"" + abs_binary_gmsh + "\" "+ abs_geo_file + " -3 -order 2 -format msh22")
        completed_process = subprocess.run([abs_binary_gmsh, abs_geo_file, "-3", "-order", "2", "-format", "msh22"], \
            check=False)
        if completed_process.returncode != 0:
            print("ERROR while running Gmsh. Mesh not created, aborting.")
            return None
        msh_filename = os.path.splitext(abs_geo_file)[0] + ".msh"
        print("Executing ", BinaryElmer, " 14 2 " + msh_filename + " -autoclean")
        completed_process = subprocess.run([BinaryElmer, "14", "2", msh_filename, "-autoclean"], \
            check=False)
        if completed_process.returncode != 0:
            print("ERROR while running ElmerGrid. Mesh not created, aborting.")
            return None
        os.remove(msh_filename)
    if do_need_meshing(abs_mesh_folder):
        print("ERROR: Failed to create mesh files. Aborting.")
        return None
    diels_file = os.path.join(directory, "diels.dat")
    shutil.copy2(diels_file, abs_mesh_folder)
    return abs_mesh_folder

# V0 is in [kV] and Vthgem1 (on the divider) in [V]
# return absolute mesh folder and absolute path to field map file (abs_mesh_folder, abs_case_result)
def solve_fields(V0, Vth1, recalculateField=False, recalculateMesh=False):
    path = os.path.abspath(__file__)
    dir_path = os.path.dirname(path)
    abs_geo_file = os.path.join(dir_path, GeoFile)
    #Detector constants:
    THGEM1_Rbot = 1.0
    THGEM1_Rtop = 1.0+8.6
    THGEM1_Rtotal = 1.0+8.6+1.2

    R1=80.0
    R2=40.0
    R3=10.0
    R4=600.0
    Rtot = R1 + 3*R2 + R3 + R4

    #distances [cm] from THGEM0 electrode to THGEM1 and cathode in the real detector
    EL_GAP_FULL = 2.20
    ANODE_THGEM1_dz=0.5
    DRIFT_L = 4.8
    LAR_EPS = 1.54 #LAr dielectric constant
    #distances from THGEM1 electrode to cathode and anode in Gmsh model #TODO: exctract from geo file
    CATHODE_dZ = 0.5
    ANODE_dZ = 0.5

    #Calculating potentials for sif file
    V_top_THGEM1 = Vth1 * THGEM1_Rtop / THGEM1_Rtotal
    V_bot_THGEM1 = Vth1 * THGEM1_Rbot / THGEM1_Rtotal
    Edrift = (V_bot_THGEM1 + 1000 * V0 * R4 / Rtot) / EL_GAP_FULL
    Einduction = V_top_THGEM1 / ANODE_THGEM1_dz
    Vcathode = V_bot_THGEM1 -  Edrift * CATHODE_dZ
    Vanode = V_top_THGEM1 - Einduction * ANODE_dZ

    abs_mesh_folder = create_mesh(abs_geo_file, recalculate=recalculateMesh)
    if abs_mesh_folder is None:
        return None
    print("Mesh folder = \"" + abs_mesh_folder + "\"")

    VERSION = str(round(Vth1)) + "v"
    abs_elmer_folder = os.path.join(dir_path, ElmerFolder)
    abs_result_file = os.path.join(abs_elmer_folder, "case_" + VERSION + ".result")
    if not os.path.isfile(abs_result_file) or recalculateField:
        abs_sif_template = os.path.join(dir_path, "case_template.sif")
        abs_sif_file = os.path.join(abs_elmer_folder, "case_" + VERSION + ".sif")

        rel_sif_file = os.path.relpath(abs_sif_file, abs_elmer_folder)
        rel_result_file = os.path.relpath(abs_result_file, abs_elmer_folder)
        os.makedirs(abs_elmer_folder, exist_ok=True)
        shutil.copy2(abs_sif_template, abs_sif_file)

        for line in fileinput.input([abs_sif_file], inplace=True):
            l = line.replace('MESH_FOLDER', abs_mesh_folder)
            l = l.replace('CASE_FILE', rel_sif_file) # Elmer does not work with absolute paths well for some reason
            l = l.replace('RESULT_FILE', rel_result_file)
            l = l.replace('CATHODE_POTENTIAL', str(Vcathode))
            l = l.replace('ANODE_POTENTIAL', str(Vanode))
            l = l.replace('TOP_THGEM_POTENTIAL', str(V_top_THGEM1))
            print(l.replace('BOT_THGEM_POTENTIAL', str(V_bot_THGEM1)), end='')

        if not os.path.isfile(abs_sif_file):
            print ("ERROR: \""+abs_sif_file+"\" file does not exist. Cannot calculate field map, aborting.")
            return None
        print("Executing \"" + BinarySolver +"\" \"" +abs_sif_file + "\"")
        print("\tfrom directory \""+abs_elmer_folder+"\"")
        completed_process = subprocess.run([BinarySolver, abs_sif_file], \
            check=False, cwd=abs_elmer_folder)
        if completed_process.returncode != 0:
            print("ERROR while running ElmerSolver. Field map is not created, aborting.")
            return None
        if not os.path.isfile(abs_result_file):
            print ("ERROR: \""+abs_result_file+"\" field map does not exist, aborting.")
            return None
    return (abs_mesh_folder, abs_result_file)


if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) < 2:
        print("Error: two to four arguments (V0, Vthgem1, force_recalculation_field, force_recalculation_mesh) must be passed.")
        sys.exit(1)
    if len(args) > 4:
        print("Warninig: Extra arguments besides (" + args[0] + " " + args[1] + " " + args[2] + " " + args[3] + ") are ignored.")
    recFeild = True
    recMesh = True
    if len(args) == 4:
        recMesh = not (args[3].lower() in ['false', '0', 'n', 'no'])
    if len(args) >= 3:
        recFeild = not (args[2].lower() in ['false', '0', 'n', 'no'])
    if len(args) >= 2:
        V0 = 0
        Vth1 = 0
        try:
            V0 = float(args[0])
        except ValueError:
            print ("Error: ", args[0], " is not a float")
            sys.exit(2)
        try:
            Vth1 = float(args[1])
        except ValueError:
            print ("Error: ", args[1], " is not a float")
            sys.exit(2)

    print("Executing \"ls\"")
    completed_process = subprocess.run(["ls"], check=False)
    if completed_process.returncode != 0:
        print("Error while running ls.")

    res = solve_fields(V0, Vth1, recalculateField=recFeild, recalculateMesh=recMesh)
    if res is None:
        sys.exit(3)
    print("Mesh folder = \"" + res[0] +"\"")
    print("Field map = \"" + res[1] + "\"")
