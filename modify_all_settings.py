#!/usr/bin/env python3
import sys
sys.dont_write_bytecode = True
import shutil
from pathlib import Path as path

""" This is a helpful script to make changes to all xml settings in this project.
This is required if some mandatory setting was added (e.g. filename with new material property)
and it is too bothersome to code it to be optional in c++ code. (The logic of when if is mandatory
and when it is not can be quite complicated.)

Another application is to replace some settings values to new ones.
It can be done only for subset of folders using filename and content filtering.
"""

def find_all_settings(root_dir: path, fname_filter = lambda fname : True, content_filter = lambda opened_file : True):
    file_list = [p for p in root_dir.rglob("*.xml") if fname_filter(p)]
    file_list_by_content = []
    for f in file_list:
        with open(f) as file:
            if content_filter(file):
                file_list_by_content.append(f)
    return file_list_by_content

def insert_text_after_marker(fname: path, markers : list, texts_to_insert : list = []):
    with open(fname, 'r') as file:
        content = file.read()
    for marker, text in zip(markers, texts_to_insert):
        content = content.replace(marker, marker + text)
    with open(fname, 'w') as file:
        file.write(content)

def replace_text(fname: path, replace_what : list, replace_with : list):
    with open(fname, 'r') as file:
        content = file.read()
    for marker, text in zip(replace_what, replace_with):
        content = content.replace(marker, text)
    with open(fname, 'w') as file:
        file.write(content)

def find_templates_only(fname: path):
    is_cache = False
    for part in fname.parts:
        if 'cache' in part:
            is_cache = True
            break
    return 'settings_' in fname.name and not is_cache

def find_no_FSQ_RG15_settings(opened_file):
    content = opened_file.read()
    return ('<FSQ_RG715_SigmaAlpha_deg>' not in content
        and '<FSQ_RG715_absorption_length_filename>' not in content
        and '<FSQ_RG715_rindex_filename>' not in content)

settings_files = find_all_settings(path('./'), fname_filter = find_templates_only, content_filter=find_no_FSQ_RG15_settings)
#print(settings_files)
markers = ["<PMMA_SigmaAlpha_deg>5</PMMA_SigmaAlpha_deg> <!-- OPT, default is 0 -->\n",
           "<TPB_emission_spectrum_filename>WLS/TPB_in_polystyrene/Emission_spectrum.dat</TPB_emission_spectrum_filename>\n"]
additions = ["    <FSQ_RG715_SigmaAlpha_deg>2</FSQ_RG715_SigmaAlpha_deg> <!-- OPT, default is 0 -->\n",
           "    <FSQ_RG715_absorption_length_filename>absorption_length/fsq_rg715_absorption_length_eV_mm_manuf.dat</FSQ_RG715_absorption_length_filename>\n" + 
           "    <FSQ_RG715_rindex_filename>refractive_index/fsq_rg715_rindex_eV_manuf.dat</FSQ_RG715_rindex_filename>\n"]
for sett in settings_files:
    insert_text_after_marker(sett, markers, texts_to_insert = additions)