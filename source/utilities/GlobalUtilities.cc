#include <utilities/GlobalUtilities.hh>

std::string int_to_str(int num)
{
  return boost::lexical_cast<std::string>(num);
}

std::string int_to_str(std::size_t num)
{
  return boost::lexical_cast<std::string>(num);
}

std::string int_to_str(int num, std::size_t decimals)
{
  std::string out = boost::lexical_cast<std::string>(num);
  if (num < 0) {
    while ((out.size()-1)<decimals) {
      out = "0" + out;
    }
  } else {
    while (out.size()<decimals) {
      out = "0" + out;
    }
  }
  return out;
}

std::string int_to_str(std::size_t num, std::size_t decimals)
{
  std::string out = boost::lexical_cast<std::string>(num);
  while (out.size()<decimals) {
    out = "0" + out;
  }
  return out;
}

std::string dbl_to_str (double val, int precision)
{
  std::stringstream ss;
  ss<<std::fixed<<std::setprecision(precision)<<val;
  return ss.str();
}

boost::optional<G4PhysicalVolumesSearchScene::Findings> FindSinglePV(std::string name)
{
  std::vector<G4PhysicalVolumesSearchScene::Findings> findingsVector = FindAllPVs(name);
  if (findingsVector.size()==0) {
    std::cout<<"********************************"<<std::endl;
    std::cout<<"FindSinglePV: could not find physical volume \""<<name<<"\"!"<<std::endl;
    return boost::none;
  }
  if (findingsVector.size()>1) {
    std::cout<<"********************************"<<std::endl;
    std::cout<<"FindSinglePV:Warning: found several physical volumes with name \""<<name<<"\"!"<<std::endl;
    std::cout<<"\treturning only the first one."<<std::endl;
  }
  return findingsVector[0];
}

std::vector<G4PhysicalVolumesSearchScene::Findings> FindAllPVs(std::string name)
{ //https://geant4-forum.web.cern.ch/t/get-volume-positioning/7403/8
  std::vector<G4PhysicalVolumesSearchScene::Findings> results;
  G4TransportationManager* transportationManager = G4TransportationManager::GetTransportationManager();
  if (nullptr == transportationManager) {
    std::cout<<"********************************"<<std::endl;
    std::cout<<"FindAllPV:Error: null G4TransportationManager."<<std::endl;
    std::cout<<"\tThis function must be called after Geant4 initialization!"<<std::endl;
    return results;
  }
  size_t nWorlds = transportationManager->GetNoWorlds();
  std::vector<G4PhysicalVolumesSearchScene::Findings> findingsVector;
  std::vector<G4VPhysicalVolume*>::iterator iterWorld = transportationManager->GetWorldsIterator();
  for (size_t i = 0; i < nWorlds; ++i, ++iterWorld) {
    G4PhysicalVolumeModel searchModel(*iterWorld); // Unlimited depth.
    G4ModelingParameters mp; // Default - no culling.
    searchModel.SetModelingParameters(&mp);
    G4PhysicalVolumesSearchScene searchScene(&searchModel, name);
    searchModel.DescribeYourselfTo(searchScene); // Initiate search.
    for (const auto& findings: searchScene.GetFindings())
      findingsVector.push_back(findings);
  }
  return findingsVector;
}

std::string strtoken(std::string &in, std::string break_symbs)
{
  std::string out;
  while (!in.empty())
  {
    char a = in.front();
    in.erase(in.begin());
    bool break_ = false;
    for (auto h = break_symbs.begin(); h != break_symbs.end(); ++h)
      if (a == *h) {
        break_ = true;
        break;
      }
    if ((break_) && (out.empty()))
      continue;
    if (break_)
      return out;
    out.push_back(a);
  }
  return out;
}

void open_output_file(std::string name, std::ofstream &str, std::ios_base::openmode _mode)
{
  ensure_file(name);
  str.open(name.c_str(), _mode);
  if (!str.is_open()){
    std::cout << "Failed to open \"" << name << "\"" << std::endl;
  }
}

void ensure_file(std::string fname)
{
  std::string folder = fname;
  while ((folder.back() != '\\') &&(folder.back()!='/') &&!folder.empty())
    folder.pop_back();
  if (!folder.empty())
    folder.pop_back();
  ensure_folder(folder);
}

void ensure_folder(std::string folder)
{
  if (folder.empty())
    return;
#if defined(_WIN32)||defined(_WIN64)
  DWORD ftyp = GetFileAttributesA(folder.c_str());
  if (!(ftyp & FILE_ATTRIBUTE_DIRECTORY) || ftyp == INVALID_FILE_ATTRIBUTES) {
    int code = system(("mkdir \"" + folder + "\"").c_str());
    if (code)
      std::cout << "mkdir error: " << GetLastError() << std::endl;
  }
#else //defined(_WIN32)||defined(_WIN64)
  struct stat st;
  if (-1==stat(folder.c_str(), &st)) {
    int err = errno;
    switch (err) {
    case (EACCES): {
      std::cout<<"Access error"<<std::endl;
      break;
    }
    case (ENAMETOOLONG): {
      std::cout<<"Path is too long"<<std::endl;
      break;
    }
    case (ENOENT) :
    case (ENOTDIR): {
      int code = system(("mkdir -p \"" + folder + "\"").c_str());
      if (code)
        std::cout << "mkdir -p error: " << code << std::endl;
      break;
    }
    default:{
      std::cout<<"stat(\""<<folder<<"\") returned -1; errno == "<<err<<std::endl;
      break;
    }
    }
  } else {
    if (!S_ISDIR(st.st_mode)) {
      int code = system(("mkdir -p \"" + folder + "\"").c_str());
      if (code)
        std::cout << "mkdir -p error: " << code << std::endl;
    }
  }
#endif //defined(_WIN32)||defined(_WIN64)
}

void rename_file(std::string origin, std::string destination)
{
#if defined(__WIN32__)
  //TODO: implement
#else
  ensure_file(destination);
  int code = system(("rm -f \"" + destination + "\"").c_str());
  code = system(("mv \"" + origin + "\" \"" + destination + "\"").c_str());
  if (code) {
    std::cerr << "mv error: " << code << std::endl;
    std::cerr << "Error renaming \""<<origin<<"\" to \""<<destination<<"\"" << std::endl;
  }
#endif //_WIN32__
}

void copy_file(std::string origin, std::string destination)
{
#if defined(__WIN32__)
  //TODO: implement
#else
  ensure_file(destination);
  int code = system(("rm -f \"" + destination + "\"").c_str());
  code = system(("cp \"" + origin + "\" " + destination).c_str());
  if (code) {
    std::cerr << "cp error: " << code << std::endl;
    std::cerr << "Error copying \""<<origin<<"\" to \""<<destination<<"\"" << std::endl;
  }
#endif //_WIN32__
}

void delete_file(std::string filename)
{
#if defined(__WIN32__)
  //TODO: implement
#else
  int code = system(("rm -f \"" + filename + "\"").c_str());
  if (code) {
    std::cerr << "rm -f error: " << code << std::endl;
    std::cerr << "Error removing \""<<filename<<"\"" << std::endl;
  }
#endif //_WIN32__
}

char* c_str_cp (const std::string &str)
{
  std::size_t i_end_= str.size();
  char* ret = new char [i_end_+1];
  for (std::size_t i=0; i!=i_end_; ++i) {
    ret[i] = str[i];
  }
  ret[i_end_] = char(0);
  return ret;
}
