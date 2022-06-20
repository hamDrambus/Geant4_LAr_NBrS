/*  File for different useful functions and classes which can
 *  be used in different classes throughout project.
 */

#ifndef GLOBAL_UTILITIES_H_
#define GLOBAL_UTILITIES_H_

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#if defined(_WIN32)||defined(_WIN64)
#define NOMINMAX
#ifndef _NO_CERN_ROOT
#include "Windows4Root.h"
#else //_NO_CERN_ROOT
#include <Windows.h>
#endif //_NO_CERN_ROOT
#include <direct.h>
#else
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#endif

#include <boost/optional.hpp>
#include <boost/lexical_cast.hpp>

#include "G4TransportationManager.hh"
#include "G4PhysicalVolumesSearchScene.hh"

boost::optional<G4PhysicalVolumesSearchScene::Findings> FindSinglePV(std::string name);
std::vector<G4PhysicalVolumesSearchScene::Findings> FindAllPVs(std::string name);

std::string int_to_str(int num);
std::string int_to_str(std::size_t num);
std::string int_to_str(int num, std::size_t decimals);
std::string int_to_str(std::size_t num, std::size_t decimals);
std::string dbl_to_str (double val, int precision=0);
std::string strtoken(std::string &in, std::string break_symbs);
void open_output_file(std::string name, std::ofstream &str, std::ios_base::openmode _mode);
void ensure_file(std::string fname); //makes sure file can be created later on
void ensure_folder(std::string folder);
void rename_file(std::string origin, std::string destination);
void copy_file(std::string origin, std::string destination);
void delete_file(std::string filename);
char* c_str_cp (const std::string &str);

#endif //GLOBAL_UTILITIES_H_
