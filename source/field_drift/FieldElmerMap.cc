#include <field_drift/FieldElmerMap.hh>

FieldElmerMap::FieldElmerMap() :
    nElements(-1),
    lastElement(-1),
    nNodes(0),
    nMaterials(0),
    hasBoundingBox(false),
    checkMultipleElement(false),
    warning(false),
    fTolerance(1e-20)
{
  elements.clear();
  nodes.clear();
  materials.clear();
  m_className = "FieldElmerMap";
  ready = false;
  debug = false;
  mapxmin = mapymin = mapzmin = 0.0;
  mapxmax = mapymax = mapzmax = 0.0;
  mapvmin = mapvmax = 0.0;
}

FieldElmerMap::FieldElmerMap(std::string header, std::string elist,
                               std::string nlist,
                               std::string volt, std::string unit) :
    nElements(-1),
    lastElement(-1),
    nNodes(-1),
    nMaterials(0),
    hasBoundingBox(false),
    checkMultipleElement(false),
    warning(false),
    fTolerance(1e-20)
{

  m_className = "FieldElmerMap";
  Initialise(header, elist, nlist, volt, unit);
}

FieldElmerMap::~FieldElmerMap()
{}

int FieldElmerMap::ReadInteger(char* token, int def, bool& error) {

  if (!token) {
    error = true;
    return def;
  }
  return atoi(token);
}

double FieldElmerMap::ReadDouble(char* token, double def, bool& error) {

  if (!token) {
    error = true;
    return def;
  }
  return atof(token);
}

void FieldElmerMap::SetMedium(const int imat, DriftMedium* m) {

  if (imat < 0 || imat >= nMaterials) {
    std::cerr << m_className << "::SetMedium:\n";
    std::cerr << "    Material index " << imat << " is out of range.\n";
    return;
  }
  if (nullptr == m) {
    if (debug) {
      std::cerr << m_className << "::SetMedium:\n";
      std::cerr << "    Medium pointer is null.\n";
    }
    materials[imat].medium = m;
    materials[imat].driftmedium = false;
    return;
  }
  if (debug) {
    std::string name = m->GetName();
    std::cout << m_className << "::SetMedium:\n";
    std::cout << "    Associated material " << imat << " with medium " << name
              << ".\n";
  }
  materials[imat].medium = m;
  materials[imat].driftmedium = m->IsDriftable();
}

DriftMedium* FieldElmerMap::GetMedium(const unsigned int& imat) const {

  if (imat >= (unsigned int)nMaterials) {
    std::cerr << m_className << "::GetMedium:\n";
    std::cerr << "    Material index " << imat << " is out of range.\n";
    return nullptr;
  }

  return materials[imat].medium;
}

bool FieldElmerMap::Initialise(std::string header, std::string elist,
                                std::string nlist, std::string volt_file,
                                std::string unit)
{
  debug = false;
  ready = false;
  nMaterials = 0;

  // Keep track of the success.
  bool ok = true;

  // Buffer for reading
  const int size = 100;
  char line[size];

  // Open the header.
  std::ifstream fheader;
  fheader.open(header.c_str(), std::ios::in);
  if (fheader.fail()) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Could not open header file " << header
              << " for reading.\n";
  }

  // Temporary variables for use in file reading
  char* token = NULL;
  bool readerror = false;
  bool readstop = false;
  int il = 0;

  // Read the header to get the number of nodes and elements.
  fheader.getline(line, size, '\n');
  token = strtok(line, " ");
  nNodes = ReadInteger(token, 0, readerror);
  token = strtok(NULL, " ");
  nElements = ReadInteger(token, 0, readerror);
  std::cout << "ComponentElmer::Initialise:\n";
  std::cout << "    Read " << nNodes << " nodes and " << nElements
            << " elements from file " << header << ".\n";
  if (readerror) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Error reading file " << header << " (line " << il
              << ").\n";
    fheader.close();
    ok = false;
    return false;
  }

  // Close the header file.
  fheader.close();

  // Open the nodes list.
  std::ifstream fnodes;
  fnodes.open(nlist.c_str(), std::ios::in);
  if (fnodes.fail()) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Could not open nodes file " << nlist << " for reading.\n";
  }

  // Check the value of the unit.
  double funit;
  if (strcmp(unit.c_str(), "mum") == 0 || strcmp(unit.c_str(), "micron") == 0 ||
      strcmp(unit.c_str(), "micrometer") == 0 || strcmp(unit.c_str(), "um") == 0) {
    funit = um;
  } else if (strcmp(unit.c_str(), "mm") == 0 ||
             strcmp(unit.c_str(), "millimeter") == 0) {
    funit = mm;
  } else if (strcmp(unit.c_str(), "cm") == 0 ||
             strcmp(unit.c_str(), "centimeter") == 0) {
    funit = cm;
  } else if (strcmp(unit.c_str(), "m") == 0 ||
             strcmp(unit.c_str(), "meter") == 0) {
    funit = m;
  } else {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Unknown length unit " << unit << ".\n";
    ok = false;
    funit = mm;
  }
  if (debug) {
    std::cout << "ComponentElmer::Initialise:\n";
    std::cout << "    Unit scaling factor = " << funit << ".\n";
  }
  // Read the nodes from the file.
  node newNode;
  newNode.w.clear();
  for (il = 0; il < nNodes; il++) {

    // Get a line from the nodes file.
    fnodes.getline(line, size, '\n');

    // Ignore the first two characters.
    token = strtok(line, " ");
    token = strtok(NULL, " ");

    // Get the node coordinates.
    token = strtok(NULL, " ");
    double xnode = ReadDouble(token, -1, readerror);
    token = strtok(NULL, " ");
    double ynode = ReadDouble(token, -1, readerror);
    token = strtok(NULL, " ");
    double znode = ReadDouble(token, -1, readerror);
    if (readerror) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Error reading file " << nlist << " (line " << il
                << ").\n";
      fnodes.close();
      ok = false;
      return false;
    }

    // Set up and create a new node.
    newNode.x = xnode * funit;
    newNode.y = ynode * funit;
    newNode.z = znode * funit;
    nodes.push_back(newNode);
  }

  // Close the nodes file.
  fnodes.close();

  // Open the potential file.
  std::ifstream fvolt;
  fvolt.open(volt_file.c_str(), std::ios::in);
  if (fvolt.fail()) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Could not open result file " << volt_file << " for reading.\n";
  }

  // Reset the line counter.
  il = 1;

  // Read past the header.
  while (!readstop && fvolt.getline(line, size, '\n')) {
    token = strtok(line, " ");
    if (strcmp(token, "Perm:") == 0) readstop = true;
    il++;
  }

  // Should have stopped: if not, print error message.
  if (!readstop) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Error reading past header of potentials file " << volt_file
              << ".\n";
    fvolt.close();
    ok = false;
    return false;
  }

  // Read past the permutation map (number of lines = nNodes).
  for (int tl = 0; tl < nNodes; tl++) {
    fvolt.getline(line, size, '\n');
    il++;
  }

  // Read the potentials.
  for (int tl = 0; tl < nNodes; tl++) {
    double v;
    fvolt.getline(line, size, '\n');
    token = strtok(line, " ");
    v = ReadDouble(token, -1, readerror);
    if (readerror) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Error reading file " << volt_file << " (line " << il
                << ").\n";
      fvolt.close();
      ok = false;
      return false;
    }
    // Place the voltage in its appropriate node. Geant4 units are used!
    nodes[tl].v = v * volt;
  }

  // Close the potentials file.
  fvolt.close();

  // Open the elements file.
  std::ifstream felems;
  felems.open(elist.c_str(), std::ios::in);
  if (felems.fail()) {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr << "    Could not open result file " << elist
              << " for reading.\n";
  }

  // Read the elements and their material indices.
  elements.clear();
  int highestnode = 0;
  element newElement;
  for (il = 0; il < nElements; il++) {

    // Get a line
    felems.getline(line, size, '\n');

    // Split into tokens.
    token = strtok(line, " ");
    // Read the 2nd-order element
    // Note: Ordering of Elmer elements can be described in the
    // ElmerSolver manual (appendix D. at the time of this comment)
    // If the order read below is compared to the shape functions used
    // eg. in ElectricField, the order is wrong, but note at the
    // end of this function the order of elements 5,6,7 will change to
    // 7,5,6 when actually recorded in newElement.emap to correct for this
    token = strtok(NULL, " ");
    int imat = ReadInteger(token, -1, readerror) - 1;
    token = strtok(NULL, " ");
    token = strtok(NULL, " ");
    int in0 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in1 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in2 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in3 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in4 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in5 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in6 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in7 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in8 = ReadInteger(token, -1, readerror);
    token = strtok(NULL, " ");
    int in9 = ReadInteger(token, -1, readerror);

    if (debug && il < 10) {
      std::cout << "    Read nodes " << in0 << ", " << in1 << ", " << in2
                << ", " << in3 << ", ... from element " << il + 1 << " of "
                << nElements << " with mat " << imat << ".\n";
    }

    // Check synchronisation.
    if (readerror) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Error reading file " << elist << " (line " << il
                << ").\n";
      felems.close();
      ok = false;
      return false;
    }

    // Check the material number and add it to material list if necessary.
    if (imat < 0) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Out-of-range material number on file " << elist
                << " (line " << il << ").\n";
      std::cerr << "    Element: " << il << ", material: " << imat << "\n";
      std::cerr << "    nodes: (" << in0 << ", " << in1 << ", " << in2 << ", "
                << in3 << ", " << in4 << ", " << in5 << ", " << in6 << ", "
                << in7 << ", " << in8 << ", " << in9 << ")\n";
      ok = false;
    }
    while(imat >= nMaterials) {
      material mat;
      mat.driftmedium = false;
      mat.medium = nullptr;
      materials.push_back(mat);
      nMaterials++;
    }

    // Check the node numbers.
    if (in0 < 1 || in1 < 1 || in2 < 1 || in3 < 1 || in4 < 1 || in5 < 1 ||
        in6 < 1 || in7 < 1 || in8 < 1 || in9 < 1) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Found a node number < 1 on file " << elist << " (line "
                << il << ").\n";
      std::cerr << "    Element: " << il << ", material: " << imat << "\n";
      std::cerr << "    nodes: (" << in0 << ", " << in1 << ", " << in2 << ", "
                << in3 << ", " << in4 << ", " << in5 << ", " << in6 << ", "
                << in7 << ", " << in8 << ", " << in9 << ")\n";
      ok = false;
    }
    if (in0 > highestnode) highestnode = in0;
    if (in1 > highestnode) highestnode = in1;
    if (in2 > highestnode) highestnode = in2;
    if (in3 > highestnode) highestnode = in3;
    if (in4 > highestnode) highestnode = in4;
    if (in5 > highestnode) highestnode = in5;
    if (in6 > highestnode) highestnode = in6;
    if (in7 > highestnode) highestnode = in7;
    if (in8 > highestnode) highestnode = in8;
    if (in9 > highestnode) highestnode = in9;

    // These elements must not be degenerate.
    if (in0 == in1 || in0 == in2 || in0 == in3 || in0 == in4 || in0 == in5 ||
        in0 == in6 || in0 == in7 || in0 == in8 || in0 == in9 || in1 == in2 ||
        in1 == in3 || in1 == in4 || in1 == in5 || in1 == in6 || in1 == in7 ||
        in1 == in8 || in1 == in9 || in2 == in3 || in2 == in4 || in2 == in5 ||
        in2 == in6 || in2 == in7 || in2 == in8 || in2 == in9 || in3 == in4 ||
        in3 == in5 || in3 == in6 || in3 == in7 || in3 == in8 || in3 == in9 ||
        in4 == in5 || in4 == in6 || in4 == in7 || in4 == in8 || in4 == in9 ||
        in5 == in6 || in5 == in7 || in5 == in8 || in5 == in9 || in6 == in7 ||
        in6 == in8 || in6 == in9 || in7 == in8 || in7 == in9 || in8 == in9) {
      std::cerr << m_className << "::Initialise:\n";
      std::cerr << "    Element " << il << " of file " << elist
                << " is degenerate,\n";
      std::cerr << "    no such elements allowed in this type of map.\n";
      ok = false;
    }

    newElement.degenerate = false;

    // Store the material reference.
    newElement.matmap = imat;

    // Node references
    newElement.emap[0] = in0 - 1;
    newElement.emap[1] = in1 - 1;
    newElement.emap[2] = in2 - 1;
    newElement.emap[3] = in3 - 1;
    newElement.emap[4] = in4 - 1;
    newElement.emap[7] = in5 - 1;
    newElement.emap[5] = in6 - 1;
    newElement.emap[6] = in7 - 1;
    newElement.emap[8] = in8 - 1;
    newElement.emap[9] = in9 - 1;
    elements.push_back(newElement);
  }

  // Close the elements file.
  felems.close();

  // Set the ready flag.
  if (ok) {
    ready = true;
  } else {
    std::cerr << m_className << "::Initialise:\n";
    std::cerr
        << "    Field map could not be read and can not be interpolated.\n";
    return false;
  }

  // Establish the ranges.
  SetRange();
  std::cout << m_className << "::Initialise:\n";
  std::cout << "    Caching the bounding box calculations of all elements.\n";
  CalculateElementBoundingBoxes();
  std::cout << m_className << "::Initialise:\n";
  std::cout << "    Caching the regions of all elements.\n";
  CalculateElementRegions();

  std::cout << m_className << "::Initialise:\n";
  std::cout << "    Finished.\n";

  return true;
}

void FieldElmerMap::SetRange()
{
  // Initial values
  mapxmin = mapymin = mapzmin = 0.;
  mapxmax = mapymax = mapzmax = 0.;
  mapvmin = mapvmax = 0.;

  // Make sure the required data is available.
  if (!ready || nNodes < 1) {
    std::cerr << m_className << "::SetRange:\n";
    std::cerr << "    Field map not yet set.\n";
    return;
  }
  if (nNodes < 1) {
    std::cerr << m_className << "::SetRange:\n";
    std::cerr << "    Number of nodes < 1.\n";
    return;
  }

  // Loop over the nodes.
  mapxmin = mapxmax = nodes[0].x;
  mapymin = mapymax = nodes[0].y;
  mapzmin = mapzmax = nodes[0].z;
  mapvmin = mapvmax = nodes[0].v;

  for (int i = 1; i < nNodes; i++) {
    if (mapxmin > nodes[i].x) mapxmin = nodes[i].x;
    if (mapxmax < nodes[i].x) mapxmax = nodes[i].x;
    if (mapymin > nodes[i].y) mapymin = nodes[i].y;
    if (mapymax < nodes[i].y) mapymax = nodes[i].y;
    if (mapzmin > nodes[i].z) mapzmin = nodes[i].z;
    if (mapzmax < nodes[i].z) mapzmax = nodes[i].z;
    if (mapvmin > nodes[i].v) mapvmin = nodes[i].v;
    if (mapvmax < nodes[i].v) mapvmax = nodes[i].v;
  }
  hasBoundingBox = true;
  // Display the range if requested.
  if (debug) PrintRange();
}

void FieldElmerMap::PrintRange() const
{
  std::cout << m_className << "::PrintRange:\n";
  std::cout << "        Dimensions of the elementary block\n";
  printf("            %15g < x < %-15g mm,\n", mapxmin / mm, mapxmax / mm);
  printf("            %15g < y < %-15g mm,\n", mapymin / mm, mapymax / mm);
  printf("            %15g < z < %-15g mm,\n", mapzmin / mm, mapzmax / mm);
  printf("            %15g < V < %-15g V.\n", mapvmin / volt, mapvmax / volt);
}

void FieldElmerMap::CalculateElementBoundingBoxes(void)
{
  // Do not proceed if not properly initialised.
  if (!ready) {
    std::cerr << m_className << "::CalculateElementBoundingBoxes:\n";
    std::cerr << "    Field map not yet initialised.\n";
    std::cerr << "    Bounding boxes of elements cannot be computed.\n";
    return;
  }
  regxmin = regymin = regzmin = DBL_MAX;
  regxmax = regymax = regzmax = -DBL_MAX;
  // Tolerance
  const double f = 0.2;
  // Calculate the bounding boxes of all elements
  for (int i = 0; i < nElements; ++i) {
    element& elem = elements[i];
    elem.xmin = std::min(
        std::min(nodes[elem.emap[0]].x, nodes[elem.emap[1]].x),
        std::min(nodes[elem.emap[2]].x, nodes[elem.emap[3]].x));
    elem.xmax = std::max(
        std::max(nodes[elem.emap[0]].x, nodes[elem.emap[1]].x),
        std::max(nodes[elem.emap[2]].x, nodes[elem.emap[3]].x));
    elem.ymin = std::min(
        std::min(nodes[elem.emap[0]].y, nodes[elem.emap[1]].y),
        std::min(nodes[elem.emap[2]].y, nodes[elem.emap[3]].y));
    elem.ymax = std::max(
        std::max(nodes[elem.emap[0]].y, nodes[elem.emap[1]].y),
        std::max(nodes[elem.emap[2]].y, nodes[elem.emap[3]].y));
    elem.zmin = std::min(
        std::min(nodes[elem.emap[0]].z, nodes[elem.emap[1]].z),
        std::min(nodes[elem.emap[2]].z, nodes[elem.emap[3]].z));
    elem.zmax = std::max(
        std::max(nodes[elem.emap[0]].z, nodes[elem.emap[1]].z),
        std::max(nodes[elem.emap[2]].z, nodes[elem.emap[3]].z));
    double tolx = f * (elem.xmax - elem.xmin), toly = f * (elem.ymax - elem.ymin), tolz = f * (elem.zmax - elem.zmin);
    elem.xmin -= tolx;
    elem.xmax += tolx;
    elem.ymin -= toly;
    elem.ymax += toly;
    elem.zmin -= tolz;
    elem.zmax += tolz;
    regxmin = std::min(regxmin, elem.xmin);
    regxmax = std::max(regxmax, elem.xmax);
    regymin = std::min(regymin, elem.ymin);
    regymax = std::max(regymax, elem.ymax);
    regzmin = std::min(regzmin, elem.zmin);
    regzmax = std::max(regzmax, elem.zmax);
  }
}

void FieldElmerMap::CalculateElementRegions(void)
{
  regions.clear();
  nRegions = 0;
  if (nElements <= 0)
    return;
  int N_split = round(pow(nElements / 1000.0, 1.0/3.0)); // along each axis. 1000 is a characteristic number of elements in each region.
  if (N_split > 16)
    N_split = 16;
  if (N_split <= 1)
    return;
  for (std::size_t rx = 0; rx != N_split; ++rx) {
    for (std::size_t ry = 0; ry != N_split; ++ry) {
      for (std::size_t rz = 0; rz != N_split; ++rz) {
        region reg;
        reg.xmin = regxmin + rx * (regxmax - regxmin) / N_split;
        reg.xmax = reg.xmin + (regxmax - regxmin) / N_split;
        reg.ymin = regymin + ry * (regymax - regymin) / N_split;
        reg.ymax = reg.ymin + (regymax - regymin) / N_split;
        reg.zmin = regzmin + rz * (regzmax - regzmin) / N_split;
        reg.zmax = reg.zmin + (regzmax - regzmin) / N_split;
        nRegions++;
        regions.push_back(reg);
      }
    }
  }
  for (int i = 0; i < nElements; ++i) {
    element& elem = elements[i];
    for (int r = 0; r < nRegions; ++r) {
      region& reg = regions[r];
      if (!(elem.xmin > reg.xmax || elem.xmax < reg.xmin
          || elem.ymin > reg.ymax || elem.ymax < reg.ymin
          || elem.zmin > reg.zmax || elem.zmax < reg.zmin))
        reg.indices.push_back(i);
    }
  }
}

void FieldElmerMap::ElectricField(const double x, const double y, const double z, double& ex,
                    double& ey, double& ez, DriftMedium* &medium, int& status)
{
  double v = 0.;
  ElectricField(x, y, z, ex, ey, ez, v, medium, status);
}

void FieldElmerMap::ElectricField(const double xin, const double yin,
                                   const double zin, double& ex, double& ey,
                                   double& ez, double& volt, DriftMedium* &medium, int& status) {
  // Copy the coordinates
  double x = xin, y = yin, z = zin;
  // Initial values
  ex = ey = ez = volt = 0.;
  status = 0;
  medium = nullptr;

  // Do not proceed if not properly initialised.
  if (!ready) {
    status = -10;
    std::cerr << m_className << "::ElectricField:\n";
    std::cerr << "    Field map not available for interpolation.\n";
    return;
  }

  if (warning) {
    std::cerr << m_className << "::ElectricField:\n";
    std::cerr << "    Warnings have been issued for this field map.\n";
  }

  // Find the element that contains this point
  double t1, t2, t3, t4, jac[4][4], det;
  int imap = FindElement13(x, y, z, t1, t2, t3, t4, jac, det);
  if (imap < 0) {
    if (debug) {
      std::cout << m_className << "::ElectricField:\n";
      std::cout << "    Point (" << x << ", " << y << ", " << z
                << " not in the mesh.\n";
    }
    status = -6;
    return;
  }

  if (debug) {
    std::cout << m_className << "::ElectricField:\n";
    std::cout << "    Global: (" << x << ", " << y << ", " << z << "),\n";
    std::cout << "    Local: (" << t1 << ", " << t2 << ", " << t3 << ", " << t4
              << " in element " << imap << "\n";
    std::cout
        << "      Node             x            y            z            V\n";
    for (int i = 0; i < 10; i++) {
      printf("      %-5d %12g %12g %12g %12g\n", elements[imap].emap[i],
             nodes[elements[imap].emap[i]].x, nodes[elements[imap].emap[i]].y,
             nodes[elements[imap].emap[i]].z, nodes[elements[imap].emap[i]].v);
    }
  }

  // Tetrahedral field
  volt = nodes[elements[imap].emap[0]].v * t1 * (2 * t1 - 1) +
         nodes[elements[imap].emap[1]].v * t2 * (2 * t2 - 1) +
         nodes[elements[imap].emap[2]].v * t3 * (2 * t3 - 1) +
         nodes[elements[imap].emap[3]].v * t4 * (2 * t4 - 1) +
         4 * nodes[elements[imap].emap[4]].v * t1 * t2 +
         4 * nodes[elements[imap].emap[5]].v * t1 * t3 +
         4 * nodes[elements[imap].emap[6]].v * t1 * t4 +
         4 * nodes[elements[imap].emap[7]].v * t2 * t3 +
         4 * nodes[elements[imap].emap[8]].v * t2 * t4 +
         4 * nodes[elements[imap].emap[9]].v * t3 * t4;
  ex = -(nodes[elements[imap].emap[0]].v * (4 * t1 - 1) * jac[0][1] +
         nodes[elements[imap].emap[1]].v * (4 * t2 - 1) * jac[1][1] +
         nodes[elements[imap].emap[2]].v * (4 * t3 - 1) * jac[2][1] +
         nodes[elements[imap].emap[3]].v * (4 * t4 - 1) * jac[3][1] +
         nodes[elements[imap].emap[4]].v *
             (4 * t2 * jac[0][1] + 4 * t1 * jac[1][1]) +
         nodes[elements[imap].emap[5]].v *
             (4 * t3 * jac[0][1] + 4 * t1 * jac[2][1]) +
         nodes[elements[imap].emap[6]].v *
             (4 * t4 * jac[0][1] + 4 * t1 * jac[3][1]) +
         nodes[elements[imap].emap[7]].v *
             (4 * t3 * jac[1][1] + 4 * t2 * jac[2][1]) +
         nodes[elements[imap].emap[8]].v *
             (4 * t4 * jac[1][1] + 4 * t2 * jac[3][1]) +
         nodes[elements[imap].emap[9]].v *
             (4 * t4 * jac[2][1] + 4 * t3 * jac[3][1])) /
       det;
  ey = -(nodes[elements[imap].emap[0]].v * (4 * t1 - 1) * jac[0][2] +
         nodes[elements[imap].emap[1]].v * (4 * t2 - 1) * jac[1][2] +
         nodes[elements[imap].emap[2]].v * (4 * t3 - 1) * jac[2][2] +
         nodes[elements[imap].emap[3]].v * (4 * t4 - 1) * jac[3][2] +
         nodes[elements[imap].emap[4]].v *
             (4 * t2 * jac[0][2] + 4 * t1 * jac[1][2]) +
         nodes[elements[imap].emap[5]].v *
             (4 * t3 * jac[0][2] + 4 * t1 * jac[2][2]) +
         nodes[elements[imap].emap[6]].v *
             (4 * t4 * jac[0][2] + 4 * t1 * jac[3][2]) +
         nodes[elements[imap].emap[7]].v *
             (4 * t3 * jac[1][2] + 4 * t2 * jac[2][2]) +
         nodes[elements[imap].emap[8]].v *
             (4 * t4 * jac[1][2] + 4 * t2 * jac[3][2]) +
         nodes[elements[imap].emap[9]].v *
             (4 * t4 * jac[2][2] + 4 * t3 * jac[3][2])) /
       det;
  ez = -(nodes[elements[imap].emap[0]].v * (4 * t1 - 1) * jac[0][3] +
         nodes[elements[imap].emap[1]].v * (4 * t2 - 1) * jac[1][3] +
         nodes[elements[imap].emap[2]].v * (4 * t3 - 1) * jac[2][3] +
         nodes[elements[imap].emap[3]].v * (4 * t4 - 1) * jac[3][3] +
         nodes[elements[imap].emap[4]].v *
             (4 * t2 * jac[0][3] + 4 * t1 * jac[1][3]) +
         nodes[elements[imap].emap[5]].v *
             (4 * t3 * jac[0][3] + 4 * t1 * jac[2][3]) +
         nodes[elements[imap].emap[6]].v *
             (4 * t4 * jac[0][3] + 4 * t1 * jac[3][3]) +
         nodes[elements[imap].emap[7]].v *
             (4 * t3 * jac[1][3] + 4 * t2 * jac[2][3]) +
         nodes[elements[imap].emap[8]].v *
             (4 * t4 * jac[1][3] + 4 * t2 * jac[3][3]) +
         nodes[elements[imap].emap[9]].v *
             (4 * t4 * jac[2][3] + 4 * t3 * jac[3][3])) /
       det;

  // Drift medium?
  medium = materials[elements[imap].matmap].medium;
  if (debug) {
    std::cout << m_className << "::ElectricField:\n";
    std::cout << "    Material " << elements[imap].matmap << ", drift flag "
              << materials[elements[imap].matmap].driftmedium << ", drift medium "
              << (medium == nullptr ? "NULL" : medium->GetName()) << "\n";
  }
  status = -5;
  if (materials[elements[imap].matmap].driftmedium) {
    status = 0;
  }
}

// Find the region for a given point.
int FieldElmerMap::FindRegion(const double x, const double y, const double z) {
  int out = -1;
  if (nRegions <= 0)
    return out;
  for (int r = 0; r < nRegions; ++r) {
    region& reg = regions[r];
    if (x <= reg.xmax && x >= reg.xmin
        && y <= reg.ymax && y >= reg.ymin
        && z <= reg.zmax && z >= reg.zmin) {
      return r;
    }
  }
  return out;
}

int FieldElmerMap::FindElement13(const double x, const double y,
                                 const double z, double& t1, double& t2,
                                 double& t3, double& t4, double jac[4][4],
                                 double& det)
{
  // Backup
  double jacbak[4][4];
  double detbak = 1.;
  double t1bak = 0., t2bak = 0., t3bak = 0., t4bak = 0.;
  int imapbak = -1;

  // Initial values.
  t1 = t2 = t3 = t4 = 0.;

  // Check previously used element
  int rc;
  int lastElem = lastElement.Get();
  if (lastElem > -1 && lastElem < nElements && !checkMultipleElement) {
    rc = Coordinates13(x, y, z, t1, t2, t3, t4, jac, det, lastElem);
    if (rc == 0 && t1 >= -fTolerance && t1 <= (1 + fTolerance) && t2 >= -fTolerance && t2 <= (1 + fTolerance) &&
        t3 >= -fTolerance && t3 <= (1 + fTolerance) && t4 >= -fTolerance && t4 <= (1 + fTolerance))
      return lastElem;
  }
  // Verify the count of volumes that contain the point.
  int nfound = 0;
  int imap = -1;

  // The loop below can be rewritten without lambda function but it will be slower by ~35%.
  auto process_found_element = [&] (int i) {
    rc = Coordinates13(x, y, z, t1, t2, t3, t4, jac, det, i);

    if (rc == 0 && t1 >= -fTolerance && t1 <= (1 + fTolerance) && t2 >= -fTolerance && t2 <= (1 + fTolerance) && t3 >= -fTolerance &&
        t3 <= (1 + fTolerance) && t4 >= -fTolerance && t4 <= (1 + fTolerance)) {
      ++nfound;
      imap = i;
      lastElement.Put(i);
      if (debug) {
        std::cout << m_className << "::FindElement13:\n";
        std::cout << "    Found matching element " << i << ".\n";
      }
      if (!checkMultipleElement) return i;
      for (int j = 0; j < 4; ++j) {
        for (int k = 0; k < 4; ++k)
          jacbak[j][k] = jac[j][k];
      }
      detbak = det;
      t1bak = t1;
      t2bak = t2;
      t3bak = t3;
      t4bak = t4;
      imapbak = imap;
      if (debug) {
        std::cout << m_className << "::FindElement13:\n";
        std::cout << "    Global = (" << x / mm << ", " << y / mm << ", " << z / mm << ")\n";
        std::cout << "    Local = (" << t1 << ", " << t2 << ", " << t3 << ", "
                  << t4 << ")\n";
        std::cout << "    Element = " << imap << "\n";
        std::cout << "                          Node             x            "
                     "y            z            V\n";
        for (int ii = 0; ii < 10; ++ii) {
          printf("                          %-5d %12g %12g %12g %12g\n",
                 elements[imap].emap[ii], nodes[elements[imap].emap[ii]].x,
                 nodes[elements[imap].emap[ii]].y,
                 nodes[elements[imap].emap[ii]].z,
                 nodes[elements[imap].emap[ii]].v);
        }
      }
    }
    return -1;
  };

  int r = FindRegion(x, y, z);
  // If r < 0 scan all elements, otherwise scan only elements in region.
  if (r < 0) {
    for (int i = 0; i != nElements; ++i) {
      element& e = elements[i];
      if (x < e.xmin || x > e.xmax ||
          y < e.ymin || y > e.ymax ||
          z < e.zmin || z > e.zmax)
        continue;
      int res = process_found_element(i);
      if (res >= 0)
        return res;
    }
  } else {
    region& reg = regions[r];
    for (int j = 0, j_end_ = reg.indices.size(); j != j_end_; ++j) {
      int i = reg.indices[j];
      element& e = elements[i];
      if (x < e.xmin || x > e.xmax ||
          y < e.ymin || y > e.ymax ||
          z < e.zmin || z > e.zmax)
        continue;
      int res = process_found_element(i);
      if (res >= 0)
        return res;
    }
  }

  // In checking mode, verify the tetrahedron/triangle count.
  if (checkMultipleElement) {
    if (nfound < 1) {
      if (debug) {
        std::cout << m_className << "::FindElement13:\n";
        std::cout << "    No element matching point (" << x << ", " << y << ", "
                  << z << ") found.\n";
      }
      lastElement.Put(-1);
      return -1;
    }
    if (nfound > 1) {
      std::cerr << m_className << "::FindElement13:\n";
      std::cerr << "    Found << " << nfound << " elements matching point ("
                << x << ", " << y << ", " << z << ").\n";
    }
    if (nfound > 0) {
      for (int j = 0; j < 4; ++j) {
        for (int k = 0; k < 4; ++k) jac[j][k] = jacbak[j][k];
      }
      det = detbak;
      t1 = t1bak;
      t2 = t2bak;
      t3 = t3bak;
      t4 = t4bak;
      imap = imapbak;
      lastElement.Put(imap);
      return imap;
    }
  }

  if (debug) {
    std::cout << m_className << "::FindElement13:\n";
    std::cout << "    No element matching point (" << x << ", " << y << ", "
              << z << ") found.\n";
  }
  return -1;
}

int FieldElmerMap::Coordinates13(double x, double y, double z, double& t1,
                                 double& t2, double& t3, double& t4,
                                 double jac[4][4], double& det, int imap) const
{
  if (debug) {
    std::cout << m_className << "::Coordinates13:\n";
    std::cout << "   Point (" << x << ", " << y << ", " << z << ")\n";
  }

  // Failure flag
  int ifail = 1;

  // Provisional values
  t1 = t2 = t3 = t4 = 0.;

  // Set tolerance parameter.
  double f = 0.5;

  // Make a first order approximation.
  int rc = Coordinates12(x, y, z, t1, t2, t3, t4, imap);
  if (rc > 0) {
    if (debug) {
      std::cout << m_className << "::Coordinates13:\n";
      std::cout << "    Failure to obtain linear estimate of isoparametric "
                   "coordinates\n";
      std::cout << "    in element " << imap << ".\n";
    }
    return ifail;
  }
  if (t1 < -f || t2 < -f || t3 < -f || t4 < -f || t1 > 1 + f || t2 > 1 + f ||
      t3 > 1 + f || t4 > 1 + f) {
    if (debug) {
      std::cout << m_className << "::Coordinates13:\n";
      std::cout << "    Linear isoparametric coordinates more than\n";
      std::cout << "    f (" << f << ") out of range in element " << imap
                << ".\n";
    }
    ifail = 0;
    return ifail;
  }

  // Start iteration.
  double td1 = t1, td2 = t2, td3 = t3, td4 = t4;
  if (debug) {
    std::cout << m_className << "::Coordinates13:\n";
    std::cout << "    Iteration starts at (t1,t2,t3,t4) = (" << td1 << ", "
              << td2 << ", " << td3 << ", " << td4 << ").\n";
  }
  // Loop
  bool converged = false;
  double diff[4], corr[4];
  for (int iter = 0; iter < 10; iter++) {
    if (debug) {
      std::cout << m_className << "::Coordinates13:\n";
      std::cout << "    Iteration " << iter << ":      (t1,t2,t3,t4) = (" << td1
                << ", " << td2 << ", " << td3 << ", " << td4 << ").\n";
    }
    // Re-compute the (x,y,z) position for this coordinate.
    double xr = nodes[elements[imap].emap[0]].x * td1 * (2 * td1 - 1) +
                nodes[elements[imap].emap[1]].x * td2 * (2 * td2 - 1) +
                nodes[elements[imap].emap[2]].x * td3 * (2 * td3 - 1) +
                nodes[elements[imap].emap[3]].x * td4 * (2 * td4 - 1) +
                nodes[elements[imap].emap[4]].x * 4 * td1 * td2 +
                nodes[elements[imap].emap[5]].x * 4 * td1 * td3 +
                nodes[elements[imap].emap[6]].x * 4 * td1 * td4 +
                nodes[elements[imap].emap[7]].x * 4 * td2 * td3 +
                nodes[elements[imap].emap[8]].x * 4 * td2 * td4 +
                nodes[elements[imap].emap[9]].x * 4 * td3 * td4;
    double yr = nodes[elements[imap].emap[0]].y * td1 * (2 * td1 - 1) +
                nodes[elements[imap].emap[1]].y * td2 * (2 * td2 - 1) +
                nodes[elements[imap].emap[2]].y * td3 * (2 * td3 - 1) +
                nodes[elements[imap].emap[3]].y * td4 * (2 * td4 - 1) +
                nodes[elements[imap].emap[4]].y * 4 * td1 * td2 +
                nodes[elements[imap].emap[5]].y * 4 * td1 * td3 +
                nodes[elements[imap].emap[6]].y * 4 * td1 * td4 +
                nodes[elements[imap].emap[7]].y * 4 * td2 * td3 +
                nodes[elements[imap].emap[8]].y * 4 * td2 * td4 +
                nodes[elements[imap].emap[9]].y * 4 * td3 * td4;
    double zr = nodes[elements[imap].emap[0]].z * td1 * (2 * td1 - 1) +
                nodes[elements[imap].emap[1]].z * td2 * (2 * td2 - 1) +
                nodes[elements[imap].emap[2]].z * td3 * (2 * td3 - 1) +
                nodes[elements[imap].emap[3]].z * td4 * (2 * td4 - 1) +
                nodes[elements[imap].emap[4]].z * 4 * td1 * td2 +
                nodes[elements[imap].emap[5]].z * 4 * td1 * td3 +
                nodes[elements[imap].emap[6]].z * 4 * td1 * td4 +
                nodes[elements[imap].emap[7]].z * 4 * td2 * td3 +
                nodes[elements[imap].emap[8]].z * 4 * td2 * td4 +
                nodes[elements[imap].emap[9]].z * 4 * td3 * td4;
    double sr = td1 + td2 + td3 + td4;

    // Compute the Jacobian.
    Jacobian13(imap, td1, td2, td3, td4, det, jac);
    // Compute the difference vector.
    diff[0] = 1 - sr;
    diff[1] = x - xr;
    diff[2] = y - yr;
    diff[3] = z - zr;

    // Update the estimate.
    for (int l = 0; l < 4; ++l) {
      corr[l] = 0;
      for (int k = 0; k < 4; ++k) {
        corr[l] += jac[l][k] * diff[k] / det;
      }
    }

    // Debugging
    if (debug) {
      std::cout << m_className << "::Coordinates13:\n";
      std::cout << "    Difference vector:  (1, x, y, z)  = (" << diff[0]
                << ", " << diff[1] << ", " << diff[2] << ", " << diff[3]
                << ").\n";
      std::cout << "    Correction vector:  (t1,t2,t3,t4) = (" << corr[0]
                << ", " << corr[1] << ", " << corr[2] << ", " << corr[3]
                << ").\n";
    }

    // Update the vector.
    td1 += corr[0];
    td2 += corr[1];
    td3 += corr[2];
    td4 += corr[3];

    // Check for convergence.
    if (fabs(corr[0]) < 1.0e-5 && fabs(corr[1]) < 1.0e-5 &&
        fabs(corr[2]) < 1.0e-5 && fabs(corr[3]) < 1.0e-5) {
      if (debug) {
        std::cout << m_className << "::Coordinates13:\n";
        std::cout << "    Convergence reached.\n";
      }
      converged = true;
      break;
    }
  }

  // No convergence reached.
  if (!converged) {
    double xmin, ymin, zmin, xmax, ymax, zmax;
    xmin = nodes[elements[imap].emap[0]].x;
    xmax = nodes[elements[imap].emap[0]].x;
    if (nodes[elements[imap].emap[1]].x < xmin) {
      xmin = nodes[elements[imap].emap[1]].x;
    }
    if (nodes[elements[imap].emap[1]].x > xmax) {
      xmax = nodes[elements[imap].emap[1]].x;
    }
    if (nodes[elements[imap].emap[2]].x < xmin) {
      xmin = nodes[elements[imap].emap[2]].x;
    }
    if (nodes[elements[imap].emap[2]].x > xmax) {
      xmax = nodes[elements[imap].emap[2]].x;
    }
    if (nodes[elements[imap].emap[3]].x < xmin) {
      xmin = nodes[elements[imap].emap[3]].x;
    }
    if (nodes[elements[imap].emap[3]].x > xmax) {
      xmax = nodes[elements[imap].emap[3]].x;
    }
    ymin = nodes[elements[imap].emap[0]].y;
    ymax = nodes[elements[imap].emap[0]].y;
    if (nodes[elements[imap].emap[1]].y < ymin) {
      ymin = nodes[elements[imap].emap[1]].y;
    }
    if (nodes[elements[imap].emap[1]].y > ymax) {
      ymax = nodes[elements[imap].emap[1]].y;
    }
    if (nodes[elements[imap].emap[2]].y < ymin) {
      ymin = nodes[elements[imap].emap[2]].y;
    }
    if (nodes[elements[imap].emap[2]].y > ymax) {
      ymax = nodes[elements[imap].emap[2]].y;
    }
    if (nodes[elements[imap].emap[3]].y < ymin) {
      ymin = nodes[elements[imap].emap[3]].y;
    }
    if (nodes[elements[imap].emap[3]].y > ymax) {
      ymax = nodes[elements[imap].emap[3]].y;
    }
    zmin = nodes[elements[imap].emap[0]].z;
    zmax = nodes[elements[imap].emap[0]].z;
    if (nodes[elements[imap].emap[1]].z < zmin) {
      zmin = nodes[elements[imap].emap[1]].z;
    }
    if (nodes[elements[imap].emap[1]].z > zmax) {
      zmax = nodes[elements[imap].emap[1]].z;
    }
    if (nodes[elements[imap].emap[2]].z < zmin) {
      zmin = nodes[elements[imap].emap[2]].z;
    }
    if (nodes[elements[imap].emap[2]].z > zmax) {
      zmax = nodes[elements[imap].emap[2]].z;
    }
    if (nodes[elements[imap].emap[3]].z < zmin) {
      zmin = nodes[elements[imap].emap[3]].z;
    }
    if (nodes[elements[imap].emap[3]].z > zmax) {
      zmax = nodes[elements[imap].emap[3]].z;
    }

    if (x >= xmin && x <= xmax && y >= ymin && y <= ymax && z >= zmin &&
        z <= zmax) {
      std::cout << m_className << "::Coordinates13:\n";
      std::cout << "    No convergence achieved "
                << "when refining internal isoparametric coordinates\n";
      std::cout << "    in element " << imap << " at position (" << x << ", "
                << y << ", " << z << ").\n";
      t1 = t2 = t3 = t4 = -1;
      return ifail;
    }
  }

  // Convergence reached.
  t1 = td1;
  t2 = td2;
  t3 = td3;
  t4 = td4;
  if (debug) {
    std::cout << m_className << "::Coordinates13:\n";
    std::cout << "    Convergence reached at (t1, t2, t3, t4) = (" << t1 << ", "
              << t2 << ", " << t3 << ", " << t4 << ").\n";
  }

  // For debugging purposes, show position.
  if (debug) {
    // Re-compute the (x,y,z) position for this coordinate.
    double xr = nodes[elements[imap].emap[0]].x * td1 * (2 * td1 - 1) +
                nodes[elements[imap].emap[1]].x * td2 * (2 * td2 - 1) +
                nodes[elements[imap].emap[2]].x * td3 * (2 * td3 - 1) +
                nodes[elements[imap].emap[3]].x * td4 * (2 * td4 - 1) +
                nodes[elements[imap].emap[4]].x * 4 * td1 * td2 +
                nodes[elements[imap].emap[5]].x * 4 * td1 * td3 +
                nodes[elements[imap].emap[6]].x * 4 * td1 * td4 +
                nodes[elements[imap].emap[7]].x * 4 * td2 * td3 +
                nodes[elements[imap].emap[8]].x * 4 * td2 * td4 +
                nodes[elements[imap].emap[9]].x * 4 * td3 * td4;
    double yr = nodes[elements[imap].emap[0]].y * td1 * (2 * td1 - 1) +
                nodes[elements[imap].emap[1]].y * td2 * (2 * td2 - 1) +
                nodes[elements[imap].emap[2]].y * td3 * (2 * td3 - 1) +
                nodes[elements[imap].emap[3]].y * td4 * (2 * td4 - 1) +
                nodes[elements[imap].emap[4]].y * 4 * td1 * td2 +
                nodes[elements[imap].emap[5]].y * 4 * td1 * td3 +
                nodes[elements[imap].emap[6]].y * 4 * td1 * td4 +
                nodes[elements[imap].emap[7]].y * 4 * td2 * td3 +
                nodes[elements[imap].emap[8]].y * 4 * td2 * td4 +
                nodes[elements[imap].emap[9]].y * 4 * td3 * td4;
    double zr = nodes[elements[imap].emap[0]].z * td1 * (2 * td1 - 1) +
                nodes[elements[imap].emap[1]].z * td2 * (2 * td2 - 1) +
                nodes[elements[imap].emap[2]].z * td3 * (2 * td3 - 1) +
                nodes[elements[imap].emap[3]].z * td4 * (2 * td4 - 1) +
                nodes[elements[imap].emap[4]].z * 4 * td1 * td2 +
                nodes[elements[imap].emap[5]].z * 4 * td1 * td3 +
                nodes[elements[imap].emap[6]].z * 4 * td1 * td4 +
                nodes[elements[imap].emap[7]].z * 4 * td2 * td3 +
                nodes[elements[imap].emap[8]].z * 4 * td2 * td4 +
                nodes[elements[imap].emap[9]].z * 4 * td3 * td4;
    double sr = td1 + td2 + td3 + td4;
    std::cout << m_className << "::Coordinates13:\n";
    std::cout << "    Position requested:     (" << x << ", " << y << ", " << z
              << ")\n";
    std::cout << "    Reconstructed:          (" << xr << ", " << yr << ", "
              << zr << ")\n";
    std::cout << "    Difference:             (" << x - xr << ", " << y - yr
              << ", " << z - zr << ")\n";
    std::cout << "    Checksum - 1:           " << sr - 1 << "\n";
  }

  // Success
  ifail = 0;
  return ifail;
}

int FieldElmerMap::Coordinates12(double x, double y, double z, double& t1,
                                     double& t2, double& t3, double& t4,
                                     int imap) const {

  if (debug) {
    std::cout << m_className << "::Coordinates12:\n";
    std::cout << "   Point (" << x << ", " << y << ", " << z << ") for element "
              << imap << "\n";
  }

  // Failure flag
  int ifail = 1;

  // Compute tetrahedral coordinates.
  t1 =
      (x - nodes[elements[imap].emap[1]].x) *
          ((nodes[elements[imap].emap[2]].y - nodes[elements[imap].emap[1]].y) *
               (nodes[elements[imap].emap[3]].z -
                nodes[elements[imap].emap[1]].z) -
           (nodes[elements[imap].emap[3]].y - nodes[elements[imap].emap[1]].y) *
               (nodes[elements[imap].emap[2]].z -
                nodes[elements[imap].emap[1]].z)) +
      (y - nodes[elements[imap].emap[1]].y) *
          ((nodes[elements[imap].emap[2]].z - nodes[elements[imap].emap[1]].z) *
               (nodes[elements[imap].emap[3]].x -
                nodes[elements[imap].emap[1]].x) -
           (nodes[elements[imap].emap[3]].z - nodes[elements[imap].emap[1]].z) *
               (nodes[elements[imap].emap[2]].x -
                nodes[elements[imap].emap[1]].x)) +
      (z - nodes[elements[imap].emap[1]].z) *
          ((nodes[elements[imap].emap[2]].x - nodes[elements[imap].emap[1]].x) *
               (nodes[elements[imap].emap[3]].y -
                nodes[elements[imap].emap[1]].y) -
           (nodes[elements[imap].emap[3]].x - nodes[elements[imap].emap[1]].x) *
               (nodes[elements[imap].emap[2]].y -
                nodes[elements[imap].emap[1]].y));
  t2 =
      (x - nodes[elements[imap].emap[2]].x) *
          ((nodes[elements[imap].emap[0]].y - nodes[elements[imap].emap[2]].y) *
               (nodes[elements[imap].emap[3]].z -
                nodes[elements[imap].emap[2]].z) -
           (nodes[elements[imap].emap[3]].y - nodes[elements[imap].emap[2]].y) *
               (nodes[elements[imap].emap[0]].z -
                nodes[elements[imap].emap[2]].z)) +
      (y - nodes[elements[imap].emap[2]].y) *
          ((nodes[elements[imap].emap[0]].z - nodes[elements[imap].emap[2]].z) *
               (nodes[elements[imap].emap[3]].x -
                nodes[elements[imap].emap[2]].x) -
           (nodes[elements[imap].emap[3]].z - nodes[elements[imap].emap[2]].z) *
               (nodes[elements[imap].emap[0]].x -
                nodes[elements[imap].emap[2]].x)) +
      (z - nodes[elements[imap].emap[2]].z) *
          ((nodes[elements[imap].emap[0]].x - nodes[elements[imap].emap[2]].x) *
               (nodes[elements[imap].emap[3]].y -
                nodes[elements[imap].emap[2]].y) -
           (nodes[elements[imap].emap[3]].x - nodes[elements[imap].emap[2]].x) *
               (nodes[elements[imap].emap[0]].y -
                nodes[elements[imap].emap[2]].y));
  t3 =
      (x - nodes[elements[imap].emap[3]].x) *
          ((nodes[elements[imap].emap[0]].y - nodes[elements[imap].emap[3]].y) *
               (nodes[elements[imap].emap[1]].z -
                nodes[elements[imap].emap[3]].z) -
           (nodes[elements[imap].emap[1]].y - nodes[elements[imap].emap[3]].y) *
               (nodes[elements[imap].emap[0]].z -
                nodes[elements[imap].emap[3]].z)) +
      (y - nodes[elements[imap].emap[3]].y) *
          ((nodes[elements[imap].emap[0]].z - nodes[elements[imap].emap[3]].z) *
               (nodes[elements[imap].emap[1]].x -
                nodes[elements[imap].emap[3]].x) -
           (nodes[elements[imap].emap[1]].z - nodes[elements[imap].emap[3]].z) *
               (nodes[elements[imap].emap[0]].x -
                nodes[elements[imap].emap[3]].x)) +
      (z - nodes[elements[imap].emap[3]].z) *
          ((nodes[elements[imap].emap[0]].x - nodes[elements[imap].emap[3]].x) *
               (nodes[elements[imap].emap[1]].y -
                nodes[elements[imap].emap[3]].y) -
           (nodes[elements[imap].emap[1]].x - nodes[elements[imap].emap[3]].x) *
               (nodes[elements[imap].emap[0]].y -
                nodes[elements[imap].emap[3]].y));
  t4 =
      (x - nodes[elements[imap].emap[0]].x) *
          ((nodes[elements[imap].emap[2]].y - nodes[elements[imap].emap[0]].y) *
               (nodes[elements[imap].emap[1]].z -
                nodes[elements[imap].emap[0]].z) -
           (nodes[elements[imap].emap[1]].y - nodes[elements[imap].emap[0]].y) *
               (nodes[elements[imap].emap[2]].z -
                nodes[elements[imap].emap[0]].z)) +
      (y - nodes[elements[imap].emap[0]].y) *
          ((nodes[elements[imap].emap[2]].z - nodes[elements[imap].emap[0]].z) *
               (nodes[elements[imap].emap[1]].x -
                nodes[elements[imap].emap[0]].x) -
           (nodes[elements[imap].emap[1]].z - nodes[elements[imap].emap[0]].z) *
               (nodes[elements[imap].emap[2]].x -
                nodes[elements[imap].emap[0]].x)) +
      (z - nodes[elements[imap].emap[0]].z) *
          ((nodes[elements[imap].emap[2]].x - nodes[elements[imap].emap[0]].x) *
               (nodes[elements[imap].emap[1]].y -
                nodes[elements[imap].emap[0]].y) -
           (nodes[elements[imap].emap[1]].x - nodes[elements[imap].emap[0]].x) *
               (nodes[elements[imap].emap[2]].y -
                nodes[elements[imap].emap[0]].y));
  t1 = t1 /
       ((nodes[elements[imap].emap[0]].x - nodes[elements[imap].emap[1]].x) *
            ((nodes[elements[imap].emap[2]].y -
              nodes[elements[imap].emap[1]].y) *
                 (nodes[elements[imap].emap[3]].z -
                  nodes[elements[imap].emap[1]].z) -
             (nodes[elements[imap].emap[3]].y -
              nodes[elements[imap].emap[1]].y) *
                 (nodes[elements[imap].emap[2]].z -
                  nodes[elements[imap].emap[1]].z)) +
        (nodes[elements[imap].emap[0]].y - nodes[elements[imap].emap[1]].y) *
            ((nodes[elements[imap].emap[2]].z -
              nodes[elements[imap].emap[1]].z) *
                 (nodes[elements[imap].emap[3]].x -
                  nodes[elements[imap].emap[1]].x) -
             (nodes[elements[imap].emap[3]].z -
              nodes[elements[imap].emap[1]].z) *
                 (nodes[elements[imap].emap[2]].x -
                  nodes[elements[imap].emap[1]].x)) +
        (nodes[elements[imap].emap[0]].z - nodes[elements[imap].emap[1]].z) *
            ((nodes[elements[imap].emap[2]].x -
              nodes[elements[imap].emap[1]].x) *
                 (nodes[elements[imap].emap[3]].y -
                  nodes[elements[imap].emap[1]].y) -
             (nodes[elements[imap].emap[3]].x -
              nodes[elements[imap].emap[1]].x) *
                 (nodes[elements[imap].emap[2]].y -
                  nodes[elements[imap].emap[1]].y)));
  t2 = t2 /
       ((nodes[elements[imap].emap[1]].x - nodes[elements[imap].emap[2]].x) *
            ((nodes[elements[imap].emap[0]].y -
              nodes[elements[imap].emap[2]].y) *
                 (nodes[elements[imap].emap[3]].z -
                  nodes[elements[imap].emap[2]].z) -
             (nodes[elements[imap].emap[3]].y -
              nodes[elements[imap].emap[2]].y) *
                 (nodes[elements[imap].emap[0]].z -
                  nodes[elements[imap].emap[2]].z)) +
        (nodes[elements[imap].emap[1]].y - nodes[elements[imap].emap[2]].y) *
            ((nodes[elements[imap].emap[0]].z -
              nodes[elements[imap].emap[2]].z) *
                 (nodes[elements[imap].emap[3]].x -
                  nodes[elements[imap].emap[2]].x) -
             (nodes[elements[imap].emap[3]].z -
              nodes[elements[imap].emap[2]].z) *
                 (nodes[elements[imap].emap[0]].x -
                  nodes[elements[imap].emap[2]].x)) +
        (nodes[elements[imap].emap[1]].z - nodes[elements[imap].emap[2]].z) *
            ((nodes[elements[imap].emap[0]].x -
              nodes[elements[imap].emap[2]].x) *
                 (nodes[elements[imap].emap[3]].y -
                  nodes[elements[imap].emap[2]].y) -
             (nodes[elements[imap].emap[3]].x -
              nodes[elements[imap].emap[2]].x) *
                 (nodes[elements[imap].emap[0]].y -
                  nodes[elements[imap].emap[2]].y)));
  t3 = t3 /
       ((nodes[elements[imap].emap[2]].x - nodes[elements[imap].emap[3]].x) *
            ((nodes[elements[imap].emap[0]].y -
              nodes[elements[imap].emap[3]].y) *
                 (nodes[elements[imap].emap[1]].z -
                  nodes[elements[imap].emap[3]].z) -
             (nodes[elements[imap].emap[1]].y -
              nodes[elements[imap].emap[3]].y) *
                 (nodes[elements[imap].emap[0]].z -
                  nodes[elements[imap].emap[3]].z)) +
        (nodes[elements[imap].emap[2]].y - nodes[elements[imap].emap[3]].y) *
            ((nodes[elements[imap].emap[0]].z -
              nodes[elements[imap].emap[3]].z) *
                 (nodes[elements[imap].emap[1]].x -
                  nodes[elements[imap].emap[3]].x) -
             (nodes[elements[imap].emap[1]].z -
              nodes[elements[imap].emap[3]].z) *
                 (nodes[elements[imap].emap[0]].x -
                  nodes[elements[imap].emap[3]].x)) +
        (nodes[elements[imap].emap[2]].z - nodes[elements[imap].emap[3]].z) *
            ((nodes[elements[imap].emap[0]].x -
              nodes[elements[imap].emap[3]].x) *
                 (nodes[elements[imap].emap[1]].y -
                  nodes[elements[imap].emap[3]].y) -
             (nodes[elements[imap].emap[1]].x -
              nodes[elements[imap].emap[3]].x) *
                 (nodes[elements[imap].emap[0]].y -
                  nodes[elements[imap].emap[3]].y)));
  t4 = t4 /
       ((nodes[elements[imap].emap[3]].x - nodes[elements[imap].emap[0]].x) *
            ((nodes[elements[imap].emap[2]].y -
              nodes[elements[imap].emap[0]].y) *
                 (nodes[elements[imap].emap[1]].z -
                  nodes[elements[imap].emap[0]].z) -
             (nodes[elements[imap].emap[1]].y -
              nodes[elements[imap].emap[0]].y) *
                 (nodes[elements[imap].emap[2]].z -
                  nodes[elements[imap].emap[0]].z)) +
        (nodes[elements[imap].emap[3]].y - nodes[elements[imap].emap[0]].y) *
            ((nodes[elements[imap].emap[2]].z -
              nodes[elements[imap].emap[0]].z) *
                 (nodes[elements[imap].emap[1]].x -
                  nodes[elements[imap].emap[0]].x) -
             (nodes[elements[imap].emap[1]].z -
              nodes[elements[imap].emap[0]].z) *
                 (nodes[elements[imap].emap[2]].x -
                  nodes[elements[imap].emap[0]].x)) +
        (nodes[elements[imap].emap[3]].z - nodes[elements[imap].emap[0]].z) *
            ((nodes[elements[imap].emap[2]].x -
              nodes[elements[imap].emap[0]].x) *
                 (nodes[elements[imap].emap[1]].y -
                  nodes[elements[imap].emap[0]].y) -
             (nodes[elements[imap].emap[1]].x -
              nodes[elements[imap].emap[0]].x) *
                 (nodes[elements[imap].emap[2]].y -
                  nodes[elements[imap].emap[0]].y)));

  // Result
  if (debug) {
    std::cout << m_className << "::Coordinates12:\n";
    std::cout << "    Tetrahedral coordinates (t, u, v, w) = (" << t1 << ", "
              << t2 << ", " << t3 << ", " << t4
              << ") sum = " << t1 + t2 + t3 + t4 << ".\n";
  }
  // Re-compute the (x,y,z) position for this coordinate.
  if (debug) {
    double xr = nodes[elements[imap].emap[0]].x * t1 +
                nodes[elements[imap].emap[1]].x * t2 +
                nodes[elements[imap].emap[2]].x * t3 +
                nodes[elements[imap].emap[3]].x * t4;
    double yr = nodes[elements[imap].emap[0]].y * t1 +
                nodes[elements[imap].emap[1]].y * t2 +
                nodes[elements[imap].emap[2]].y * t3 +
                nodes[elements[imap].emap[3]].y * t4;
    double zr = nodes[elements[imap].emap[0]].z * t1 +
                nodes[elements[imap].emap[1]].z * t2 +
                nodes[elements[imap].emap[2]].z * t3 +
                nodes[elements[imap].emap[3]].z * t4;
    double sr = t1 + t2 + t3 + t4;
    std::cout << m_className << "::Coordinates12:\n";
    std::cout << "    Position requested:     (" << x << ", " << y << ", " << z
              << ")\n";
    std::cout << "    Reconstructed:          (" << xr << ", " << yr << ", "
              << zr << ")\n";
    std::cout << "    Difference:             (" << x - xr << ", " << y - yr
              << ", " << z - zr << ")\n";
    std::cout << "    Checksum - 1:           " << sr - 1 << "\n";
  }

  // This should always work.
  ifail = 0;
  return ifail;
}

void FieldElmerMap::Jacobian13(int i, double t, double u, double v,
                               double w, double& det, double jac[4][4]) const {

  // Initial values
  det = 0;
  for (int j = 0; j < 4; ++j) {
    for (int k = 0; k < 4; ++k) jac[j][k] = 0;
  }

  // Be sure that the element is within range
  if (i < 0 || i >= nElements) {
    std::cerr << m_className << "::Jacobian13:\n";
    std::cerr << "    Element " << i << " out of range.\n";
    return;
  }

  // Determinant of the quadrilateral serendipity Jacobian
  det =
      -(((-4 * v * nodes[elements[i].emap[9]].x - nodes[elements[i].emap[1]].x +
          4 * u * nodes[elements[i].emap[1]].x + nodes[elements[i].emap[3]].x -
          4 * w * nodes[elements[i].emap[3]].x +
          4 * t * nodes[elements[i].emap[4]].x -
          4 * t * nodes[elements[i].emap[6]].x +
          4 * v * nodes[elements[i].emap[7]].x -
          4 * u * nodes[elements[i].emap[8]].x +
          4 * w * nodes[elements[i].emap[8]].x) *
             (4 * w * nodes[elements[i].emap[9]].y -
              nodes[elements[i].emap[2]].y +
              4 * v * nodes[elements[i].emap[2]].y +
              4 * t * nodes[elements[i].emap[5]].y +
              4 * u * nodes[elements[i].emap[7]].y) +
         (nodes[elements[i].emap[1]].x - 4 * u * nodes[elements[i].emap[1]].x -
          nodes[elements[i].emap[2]].x + 4 * v * nodes[elements[i].emap[2]].x -
          4 * t * nodes[elements[i].emap[4]].x +
          4 * t * nodes[elements[i].emap[5]].x +
          4 * u * nodes[elements[i].emap[7]].x -
          4 * v * nodes[elements[i].emap[7]].x +
          4 * w *
              (nodes[elements[i].emap[9]].x - nodes[elements[i].emap[8]].x)) *
             (4 * v * nodes[elements[i].emap[9]].y -
              nodes[elements[i].emap[3]].y +
              4 * w * nodes[elements[i].emap[3]].y +
              4 * t * nodes[elements[i].emap[6]].y +
              4 * u * nodes[elements[i].emap[8]].y) +
         (-4 * w * nodes[elements[i].emap[9]].x +
          4 * v *
              (nodes[elements[i].emap[9]].x - nodes[elements[i].emap[2]].x) +
          nodes[elements[i].emap[2]].x - nodes[elements[i].emap[3]].x +
          4 * w * nodes[elements[i].emap[3]].x -
          4 * t * nodes[elements[i].emap[5]].x +
          4 * t * nodes[elements[i].emap[6]].x -
          4 * u * nodes[elements[i].emap[7]].x +
          4 * u * nodes[elements[i].emap[8]].x) *
             ((-1 + 4 * u) * nodes[elements[i].emap[1]].y +
              4 * (t * nodes[elements[i].emap[4]].y +
                   v * nodes[elements[i].emap[7]].y +
                   w * nodes[elements[i].emap[8]].y))) *
        ((-1 + 4 * t) * nodes[elements[i].emap[0]].z +
         4 * (u * nodes[elements[i].emap[4]].z +
              v * nodes[elements[i].emap[5]].z +
              w * nodes[elements[i].emap[6]].z))) -
      ((nodes[elements[i].emap[1]].x - 4 * u * nodes[elements[i].emap[1]].x -
        nodes[elements[i].emap[3]].x + 4 * w * nodes[elements[i].emap[3]].x -
        4 * t * nodes[elements[i].emap[4]].x +
        4 * t * nodes[elements[i].emap[6]].x +
        4 * v * (nodes[elements[i].emap[9]].x - nodes[elements[i].emap[7]].x) +
        4 * u * nodes[elements[i].emap[8]].x -
        4 * w * nodes[elements[i].emap[8]].x) *
           ((-1 + 4 * t) * nodes[elements[i].emap[0]].y +
            4 * (u * nodes[elements[i].emap[4]].y +
                 v * nodes[elements[i].emap[5]].y +
                 w * nodes[elements[i].emap[6]].y)) -
       ((-1 + 4 * t) * nodes[elements[i].emap[0]].x +
        nodes[elements[i].emap[1]].x - 4 * u * nodes[elements[i].emap[1]].x +
        4 * (-(t * nodes[elements[i].emap[4]].x) +
             u * nodes[elements[i].emap[4]].x +
             v * nodes[elements[i].emap[5]].x +
             w * nodes[elements[i].emap[6]].x -
             v * nodes[elements[i].emap[7]].x -
             w * nodes[elements[i].emap[8]].x)) *
           (4 * v * nodes[elements[i].emap[9]].y -
            nodes[elements[i].emap[3]].y +
            4 * w * nodes[elements[i].emap[3]].y +
            4 * t * nodes[elements[i].emap[6]].y +
            4 * u * nodes[elements[i].emap[8]].y) +
       ((-1 + 4 * t) * nodes[elements[i].emap[0]].x -
        4 * v * nodes[elements[i].emap[9]].x + nodes[elements[i].emap[3]].x -
        4 * w * nodes[elements[i].emap[3]].x +
        4 * u * nodes[elements[i].emap[4]].x +
        4 * v * nodes[elements[i].emap[5]].x -
        4 * t * nodes[elements[i].emap[6]].x +
        4 * w * nodes[elements[i].emap[6]].x -
        4 * u * nodes[elements[i].emap[8]].x) *
           ((-1 + 4 * u) * nodes[elements[i].emap[1]].y +
            4 * (t * nodes[elements[i].emap[4]].y +
                 v * nodes[elements[i].emap[7]].y +
                 w * nodes[elements[i].emap[8]].y))) *
          (4 * w * nodes[elements[i].emap[9]].z - nodes[elements[i].emap[2]].z +
           4 * v * nodes[elements[i].emap[2]].z +
           4 * t * nodes[elements[i].emap[5]].z +
           4 * u * nodes[elements[i].emap[7]].z) +
      ((nodes[elements[i].emap[1]].x - 4 * u * nodes[elements[i].emap[1]].x -
        nodes[elements[i].emap[2]].x + 4 * v * nodes[elements[i].emap[2]].x -
        4 * t * nodes[elements[i].emap[4]].x +
        4 * t * nodes[elements[i].emap[5]].x +
        4 * u * nodes[elements[i].emap[7]].x -
        4 * v * nodes[elements[i].emap[7]].x +
        4 * w * (nodes[elements[i].emap[9]].x - nodes[elements[i].emap[8]].x)) *
           ((-1 + 4 * t) * nodes[elements[i].emap[0]].y +
            4 * (u * nodes[elements[i].emap[4]].y +
                 v * nodes[elements[i].emap[5]].y +
                 w * nodes[elements[i].emap[6]].y)) -
       ((-1 + 4 * t) * nodes[elements[i].emap[0]].x +
        nodes[elements[i].emap[1]].x - 4 * u * nodes[elements[i].emap[1]].x +
        4 * (-(t * nodes[elements[i].emap[4]].x) +
             u * nodes[elements[i].emap[4]].x +
             v * nodes[elements[i].emap[5]].x +
             w * nodes[elements[i].emap[6]].x -
             v * nodes[elements[i].emap[7]].x -
             w * nodes[elements[i].emap[8]].x)) *
           (4 * w * nodes[elements[i].emap[9]].y -
            nodes[elements[i].emap[2]].y +
            4 * v * nodes[elements[i].emap[2]].y +
            4 * t * nodes[elements[i].emap[5]].y +
            4 * u * nodes[elements[i].emap[7]].y) +
       ((-1 + 4 * t) * nodes[elements[i].emap[0]].x -
        4 * w * nodes[elements[i].emap[9]].x + nodes[elements[i].emap[2]].x -
        4 * v * nodes[elements[i].emap[2]].x +
        4 * u * nodes[elements[i].emap[4]].x -
        4 * t * nodes[elements[i].emap[5]].x +
        4 * v * nodes[elements[i].emap[5]].x +
        4 * w * nodes[elements[i].emap[6]].x -
        4 * u * nodes[elements[i].emap[7]].x) *
           ((-1 + 4 * u) * nodes[elements[i].emap[1]].y +
            4 * (t * nodes[elements[i].emap[4]].y +
                 v * nodes[elements[i].emap[7]].y +
                 w * nodes[elements[i].emap[8]].y))) *
          (4 * v * nodes[elements[i].emap[9]].z - nodes[elements[i].emap[3]].z +
           4 * w * nodes[elements[i].emap[3]].z +
           4 * t * nodes[elements[i].emap[6]].z +
           4 * u * nodes[elements[i].emap[8]].z) +
      ((-4 * w * nodes[elements[i].emap[9]].x +
        4 * v * (nodes[elements[i].emap[9]].x - nodes[elements[i].emap[2]].x) +
        nodes[elements[i].emap[2]].x - nodes[elements[i].emap[3]].x +
        4 * w * nodes[elements[i].emap[3]].x -
        4 * t * nodes[elements[i].emap[5]].x +
        4 * t * nodes[elements[i].emap[6]].x -
        4 * u * nodes[elements[i].emap[7]].x +
        4 * u * nodes[elements[i].emap[8]].x) *
           ((-1 + 4 * t) * nodes[elements[i].emap[0]].y +
            4 * (u * nodes[elements[i].emap[4]].y +
                 v * nodes[elements[i].emap[5]].y +
                 w * nodes[elements[i].emap[6]].y)) +
       ((-1 + 4 * t) * nodes[elements[i].emap[0]].x -
        4 * v * nodes[elements[i].emap[9]].x + nodes[elements[i].emap[3]].x -
        4 * w * nodes[elements[i].emap[3]].x +
        4 * u * nodes[elements[i].emap[4]].x +
        4 * v * nodes[elements[i].emap[5]].x -
        4 * t * nodes[elements[i].emap[6]].x +
        4 * w * nodes[elements[i].emap[6]].x -
        4 * u * nodes[elements[i].emap[8]].x) *
           (4 * w * nodes[elements[i].emap[9]].y -
            nodes[elements[i].emap[2]].y +
            4 * v * nodes[elements[i].emap[2]].y +
            4 * t * nodes[elements[i].emap[5]].y +
            4 * u * nodes[elements[i].emap[7]].y) -
       ((-1 + 4 * t) * nodes[elements[i].emap[0]].x -
        4 * w * nodes[elements[i].emap[9]].x + nodes[elements[i].emap[2]].x -
        4 * v * nodes[elements[i].emap[2]].x +
        4 * u * nodes[elements[i].emap[4]].x -
        4 * t * nodes[elements[i].emap[5]].x +
        4 * v * nodes[elements[i].emap[5]].x +
        4 * w * nodes[elements[i].emap[6]].x -
        4 * u * nodes[elements[i].emap[7]].x) *
           (4 * v * nodes[elements[i].emap[9]].y -
            nodes[elements[i].emap[3]].y +
            4 * w * nodes[elements[i].emap[3]].y +
            4 * t * nodes[elements[i].emap[6]].y +
            4 * u * nodes[elements[i].emap[8]].y)) *
          ((-1 + 4 * u) * nodes[elements[i].emap[1]].z +
           4 * (t * nodes[elements[i].emap[4]].z +
                v * nodes[elements[i].emap[7]].z +
                w * nodes[elements[i].emap[8]].z));

  jac[0][0] =
      -((((-1 + 4 * u) * nodes[elements[i].emap[1]].x +
          4 * (t * nodes[elements[i].emap[4]].x +
               v * nodes[elements[i].emap[7]].x +
               w * nodes[elements[i].emap[8]].x)) *
             (4 * v * nodes[elements[i].emap[9]].y -
              nodes[elements[i].emap[3]].y +
              4 * w * nodes[elements[i].emap[3]].y +
              4 * t * nodes[elements[i].emap[6]].y +
              4 * u * nodes[elements[i].emap[8]].y) -
         (4 * v * nodes[elements[i].emap[9]].x - nodes[elements[i].emap[3]].x +
          4 * w * nodes[elements[i].emap[3]].x +
          4 * t * nodes[elements[i].emap[6]].x +
          4 * u * nodes[elements[i].emap[8]].x) *
             ((-1 + 4 * u) * nodes[elements[i].emap[1]].y +
              4 * (t * nodes[elements[i].emap[4]].y +
                   v * nodes[elements[i].emap[7]].y +
                   w * nodes[elements[i].emap[8]].y))) *
        (4 * w * nodes[elements[i].emap[9]].z - nodes[elements[i].emap[2]].z +
         4 * v * nodes[elements[i].emap[2]].z +
         4 * t * nodes[elements[i].emap[5]].z +
         4 * u * nodes[elements[i].emap[7]].z)) +
      (((-1 + 4 * u) * nodes[elements[i].emap[1]].x +
        4 * (t * nodes[elements[i].emap[4]].x +
             v * nodes[elements[i].emap[7]].x +
             w * nodes[elements[i].emap[8]].x)) *
           (4 * w * nodes[elements[i].emap[9]].y -
            nodes[elements[i].emap[2]].y +
            4 * v * nodes[elements[i].emap[2]].y +
            4 * t * nodes[elements[i].emap[5]].y +
            4 * u * nodes[elements[i].emap[7]].y) -
       (4 * w * nodes[elements[i].emap[9]].x - nodes[elements[i].emap[2]].x +
        4 * v * nodes[elements[i].emap[2]].x +
        4 * t * nodes[elements[i].emap[5]].x +
        4 * u * nodes[elements[i].emap[7]].x) *
           ((-1 + 4 * u) * nodes[elements[i].emap[1]].y +
            4 * (t * nodes[elements[i].emap[4]].y +
                 v * nodes[elements[i].emap[7]].y +
                 w * nodes[elements[i].emap[8]].y))) *
          (4 * v * nodes[elements[i].emap[9]].z - nodes[elements[i].emap[3]].z +
           4 * w * nodes[elements[i].emap[3]].z +
           4 * t * nodes[elements[i].emap[6]].z +
           4 * u * nodes[elements[i].emap[8]].z) +
      (-((4 * v * nodes[elements[i].emap[9]].x - nodes[elements[i].emap[3]].x +
          4 * w * nodes[elements[i].emap[3]].x +
          4 * t * nodes[elements[i].emap[6]].x +
          4 * u * nodes[elements[i].emap[8]].x) *
         (4 * w * nodes[elements[i].emap[9]].y - nodes[elements[i].emap[2]].y +
          4 * v * nodes[elements[i].emap[2]].y +
          4 * t * nodes[elements[i].emap[5]].y +
          4 * u * nodes[elements[i].emap[7]].y)) +
       (4 * w * nodes[elements[i].emap[9]].x - nodes[elements[i].emap[2]].x +
        4 * v * nodes[elements[i].emap[2]].x +
        4 * t * nodes[elements[i].emap[5]].x +
        4 * u * nodes[elements[i].emap[7]].x) *
           (4 * v * nodes[elements[i].emap[9]].y -
            nodes[elements[i].emap[3]].y +
            4 * w * nodes[elements[i].emap[3]].y +
            4 * t * nodes[elements[i].emap[6]].y +
            4 * u * nodes[elements[i].emap[8]].y)) *
          ((-1 + 4 * u) * nodes[elements[i].emap[1]].z +
           4 * (t * nodes[elements[i].emap[4]].z +
                v * nodes[elements[i].emap[7]].z +
                w * nodes[elements[i].emap[8]].z));

  jac[0][1] =
      (nodes[elements[i].emap[1]].y - 4 * u * nodes[elements[i].emap[1]].y -
       nodes[elements[i].emap[3]].y + 4 * w * nodes[elements[i].emap[3]].y -
       4 * t * nodes[elements[i].emap[4]].y +
       4 * t * nodes[elements[i].emap[6]].y +
       4 * v * (nodes[elements[i].emap[9]].y - nodes[elements[i].emap[7]].y) +
       4 * u * nodes[elements[i].emap[8]].y -
       4 * w * nodes[elements[i].emap[8]].y) *
          (4 * w * nodes[elements[i].emap[9]].z - nodes[elements[i].emap[2]].z +
           4 * v * nodes[elements[i].emap[2]].z +
           4 * t * nodes[elements[i].emap[5]].z +
           4 * u * nodes[elements[i].emap[7]].z) +
      (-4 * w * nodes[elements[i].emap[9]].y - nodes[elements[i].emap[1]].y +
       4 * u * nodes[elements[i].emap[1]].y + nodes[elements[i].emap[2]].y -
       4 * v * nodes[elements[i].emap[2]].y +
       4 * t * nodes[elements[i].emap[4]].y -
       4 * t * nodes[elements[i].emap[5]].y -
       4 * u * nodes[elements[i].emap[7]].y +
       4 * v * nodes[elements[i].emap[7]].y +
       4 * w * nodes[elements[i].emap[8]].y) *
          (4 * v * nodes[elements[i].emap[9]].z - nodes[elements[i].emap[3]].z +
           4 * w * nodes[elements[i].emap[3]].z +
           4 * t * nodes[elements[i].emap[6]].z +
           4 * u * nodes[elements[i].emap[8]].z) +
      (-4 * v * nodes[elements[i].emap[9]].y +
       4 * w * nodes[elements[i].emap[9]].y - nodes[elements[i].emap[2]].y +
       4 * v * nodes[elements[i].emap[2]].y + nodes[elements[i].emap[3]].y -
       4 * w * nodes[elements[i].emap[3]].y +
       4 * t * nodes[elements[i].emap[5]].y -
       4 * t * nodes[elements[i].emap[6]].y +
       4 * u * nodes[elements[i].emap[7]].y -
       4 * u * nodes[elements[i].emap[8]].y) *
          ((-1 + 4 * u) * nodes[elements[i].emap[1]].z +
           4 * (t * nodes[elements[i].emap[4]].z +
                v * nodes[elements[i].emap[7]].z +
                w * nodes[elements[i].emap[8]].z));

  jac[0][2] =
      (-4 * v * nodes[elements[i].emap[9]].x - nodes[elements[i].emap[1]].x +
       4 * u * nodes[elements[i].emap[1]].x + nodes[elements[i].emap[3]].x -
       4 * w * nodes[elements[i].emap[3]].x +
       4 * t * nodes[elements[i].emap[4]].x -
       4 * t * nodes[elements[i].emap[6]].x +
       4 * v * nodes[elements[i].emap[7]].x -
       4 * u * nodes[elements[i].emap[8]].x +
       4 * w * nodes[elements[i].emap[8]].x) *
          (4 * w * nodes[elements[i].emap[9]].z - nodes[elements[i].emap[2]].z +
           4 * v * nodes[elements[i].emap[2]].z +
           4 * t * nodes[elements[i].emap[5]].z +
           4 * u * nodes[elements[i].emap[7]].z) +
      (nodes[elements[i].emap[1]].x - 4 * u * nodes[elements[i].emap[1]].x -
       nodes[elements[i].emap[2]].x + 4 * v * nodes[elements[i].emap[2]].x -
       4 * t * nodes[elements[i].emap[4]].x +
       4 * t * nodes[elements[i].emap[5]].x +
       4 * u * nodes[elements[i].emap[7]].x -
       4 * v * nodes[elements[i].emap[7]].x +
       4 * w * (nodes[elements[i].emap[9]].x - nodes[elements[i].emap[8]].x)) *
          (4 * v * nodes[elements[i].emap[9]].z - nodes[elements[i].emap[3]].z +
           4 * w * nodes[elements[i].emap[3]].z +
           4 * t * nodes[elements[i].emap[6]].z +
           4 * u * nodes[elements[i].emap[8]].z) +
      (-4 * w * nodes[elements[i].emap[9]].x +
       4 * v * (nodes[elements[i].emap[9]].x - nodes[elements[i].emap[2]].x) +
       nodes[elements[i].emap[2]].x - nodes[elements[i].emap[3]].x +
       4 * w * nodes[elements[i].emap[3]].x -
       4 * t * nodes[elements[i].emap[5]].x +
       4 * t * nodes[elements[i].emap[6]].x -
       4 * u * nodes[elements[i].emap[7]].x +
       4 * u * nodes[elements[i].emap[8]].x) *
          ((-1 + 4 * u) * nodes[elements[i].emap[1]].z +
           4 * (t * nodes[elements[i].emap[4]].z +
                v * nodes[elements[i].emap[7]].z +
                w * nodes[elements[i].emap[8]].z));

  jac[0][3] =
      (nodes[elements[i].emap[1]].x - 4 * u * nodes[elements[i].emap[1]].x -
       nodes[elements[i].emap[3]].x + 4 * w * nodes[elements[i].emap[3]].x -
       4 * t * nodes[elements[i].emap[4]].x +
       4 * t * nodes[elements[i].emap[6]].x +
       4 * v * (nodes[elements[i].emap[9]].x - nodes[elements[i].emap[7]].x) +
       4 * u * nodes[elements[i].emap[8]].x -
       4 * w * nodes[elements[i].emap[8]].x) *
          (4 * w * nodes[elements[i].emap[9]].y - nodes[elements[i].emap[2]].y +
           4 * v * nodes[elements[i].emap[2]].y +
           4 * t * nodes[elements[i].emap[5]].y +
           4 * u * nodes[elements[i].emap[7]].y) +
      (-4 * w * nodes[elements[i].emap[9]].x - nodes[elements[i].emap[1]].x +
       4 * u * nodes[elements[i].emap[1]].x + nodes[elements[i].emap[2]].x -
       4 * v * nodes[elements[i].emap[2]].x +
       4 * t * nodes[elements[i].emap[4]].x -
       4 * t * nodes[elements[i].emap[5]].x -
       4 * u * nodes[elements[i].emap[7]].x +
       4 * v * nodes[elements[i].emap[7]].x +
       4 * w * nodes[elements[i].emap[8]].x) *
          (4 * v * nodes[elements[i].emap[9]].y - nodes[elements[i].emap[3]].y +
           4 * w * nodes[elements[i].emap[3]].y +
           4 * t * nodes[elements[i].emap[6]].y +
           4 * u * nodes[elements[i].emap[8]].y) +
      (-4 * v * nodes[elements[i].emap[9]].x +
       4 * w * nodes[elements[i].emap[9]].x - nodes[elements[i].emap[2]].x +
       4 * v * nodes[elements[i].emap[2]].x + nodes[elements[i].emap[3]].x -
       4 * w * nodes[elements[i].emap[3]].x +
       4 * t * nodes[elements[i].emap[5]].x -
       4 * t * nodes[elements[i].emap[6]].x +
       4 * u * nodes[elements[i].emap[7]].x -
       4 * u * nodes[elements[i].emap[8]].x) *
          ((-1 + 4 * u) * nodes[elements[i].emap[1]].y +
           4 * (t * nodes[elements[i].emap[4]].y +
                v * nodes[elements[i].emap[7]].y +
                w * nodes[elements[i].emap[8]].y));

  jac[1][0] =
      -((-((4 * v * nodes[elements[i].emap[9]].x -
            nodes[elements[i].emap[3]].x +
            4 * w * nodes[elements[i].emap[3]].x +
            4 * t * nodes[elements[i].emap[6]].x +
            4 * u * nodes[elements[i].emap[8]].x) *
           (4 * w * nodes[elements[i].emap[9]].y -
            nodes[elements[i].emap[2]].y +
            4 * v * nodes[elements[i].emap[2]].y +
            4 * t * nodes[elements[i].emap[5]].y +
            4 * u * nodes[elements[i].emap[7]].y)) +
         (4 * w * nodes[elements[i].emap[9]].x - nodes[elements[i].emap[2]].x +
          4 * v * nodes[elements[i].emap[2]].x +
          4 * t * nodes[elements[i].emap[5]].x +
          4 * u * nodes[elements[i].emap[7]].x) *
             (4 * v * nodes[elements[i].emap[9]].y -
              nodes[elements[i].emap[3]].y +
              4 * w * nodes[elements[i].emap[3]].y +
              4 * t * nodes[elements[i].emap[6]].y +
              4 * u * nodes[elements[i].emap[8]].y)) *
        ((-1 + 4 * t) * nodes[elements[i].emap[0]].z +
         4 * (u * nodes[elements[i].emap[4]].z +
              v * nodes[elements[i].emap[5]].z +
              w * nodes[elements[i].emap[6]].z))) +
      (-((4 * v * nodes[elements[i].emap[9]].x - nodes[elements[i].emap[3]].x +
          4 * w * nodes[elements[i].emap[3]].x +
          4 * t * nodes[elements[i].emap[6]].x +
          4 * u * nodes[elements[i].emap[8]].x) *
         ((-1 + 4 * t) * nodes[elements[i].emap[0]].y +
          4 * (u * nodes[elements[i].emap[4]].y +
               v * nodes[elements[i].emap[5]].y +
               w * nodes[elements[i].emap[6]].y))) +
       ((-1 + 4 * t) * nodes[elements[i].emap[0]].x +
        4 * (u * nodes[elements[i].emap[4]].x +
             v * nodes[elements[i].emap[5]].x +
             w * nodes[elements[i].emap[6]].x)) *
           (4 * v * nodes[elements[i].emap[9]].y -
            nodes[elements[i].emap[3]].y +
            4 * w * nodes[elements[i].emap[3]].y +
            4 * t * nodes[elements[i].emap[6]].y +
            4 * u * nodes[elements[i].emap[8]].y)) *
          (4 * w * nodes[elements[i].emap[9]].z - nodes[elements[i].emap[2]].z +
           4 * v * nodes[elements[i].emap[2]].z +
           4 * t * nodes[elements[i].emap[5]].z +
           4 * u * nodes[elements[i].emap[7]].z) -
      (-((4 * w * nodes[elements[i].emap[9]].x - nodes[elements[i].emap[2]].x +
          4 * v * nodes[elements[i].emap[2]].x +
          4 * t * nodes[elements[i].emap[5]].x +
          4 * u * nodes[elements[i].emap[7]].x) *
         ((-1 + 4 * t) * nodes[elements[i].emap[0]].y +
          4 * (u * nodes[elements[i].emap[4]].y +
               v * nodes[elements[i].emap[5]].y +
               w * nodes[elements[i].emap[6]].y))) +
       ((-1 + 4 * t) * nodes[elements[i].emap[0]].x +
        4 * (u * nodes[elements[i].emap[4]].x +
             v * nodes[elements[i].emap[5]].x +
             w * nodes[elements[i].emap[6]].x)) *
           (4 * w * nodes[elements[i].emap[9]].y -
            nodes[elements[i].emap[2]].y +
            4 * v * nodes[elements[i].emap[2]].y +
            4 * t * nodes[elements[i].emap[5]].y +
            4 * u * nodes[elements[i].emap[7]].y)) *
          (4 * v * nodes[elements[i].emap[9]].z - nodes[elements[i].emap[3]].z +
           4 * w * nodes[elements[i].emap[3]].z +
           4 * t * nodes[elements[i].emap[6]].z +
           4 * u * nodes[elements[i].emap[8]].z);

  jac[1][1] =
      (-4 * w * nodes[elements[i].emap[9]].y +
       4 * v * (nodes[elements[i].emap[9]].y - nodes[elements[i].emap[2]].y) +
       nodes[elements[i].emap[2]].y - nodes[elements[i].emap[3]].y +
       4 * w * nodes[elements[i].emap[3]].y -
       4 * t * nodes[elements[i].emap[5]].y +
       4 * t * nodes[elements[i].emap[6]].y -
       4 * u * nodes[elements[i].emap[7]].y +
       4 * u * nodes[elements[i].emap[8]].y) *
          ((-1 + 4 * t) * nodes[elements[i].emap[0]].z +
           4 * (u * nodes[elements[i].emap[4]].z +
                v * nodes[elements[i].emap[5]].z +
                w * nodes[elements[i].emap[6]].z)) +
      ((-1 + 4 * t) * nodes[elements[i].emap[0]].y -
       4 * v * nodes[elements[i].emap[9]].y + nodes[elements[i].emap[3]].y -
       4 * w * nodes[elements[i].emap[3]].y +
       4 * u * nodes[elements[i].emap[4]].y +
       4 * v * nodes[elements[i].emap[5]].y -
       4 * t * nodes[elements[i].emap[6]].y +
       4 * w * nodes[elements[i].emap[6]].y -
       4 * u * nodes[elements[i].emap[8]].y) *
          (4 * w * nodes[elements[i].emap[9]].z - nodes[elements[i].emap[2]].z +
           4 * v * nodes[elements[i].emap[2]].z +
           4 * t * nodes[elements[i].emap[5]].z +
           4 * u * nodes[elements[i].emap[7]].z) -
      ((-1 + 4 * t) * nodes[elements[i].emap[0]].y -
       4 * w * nodes[elements[i].emap[9]].y + nodes[elements[i].emap[2]].y -
       4 * v * nodes[elements[i].emap[2]].y +
       4 * u * nodes[elements[i].emap[4]].y -
       4 * t * nodes[elements[i].emap[5]].y +
       4 * v * nodes[elements[i].emap[5]].y +
       4 * w * nodes[elements[i].emap[6]].y -
       4 * u * nodes[elements[i].emap[7]].y) *
          (4 * v * nodes[elements[i].emap[9]].z - nodes[elements[i].emap[3]].z +
           4 * w * nodes[elements[i].emap[3]].z +
           4 * t * nodes[elements[i].emap[6]].z +
           4 * u * nodes[elements[i].emap[8]].z);

  jac[1][2] =
      (-4 * v * nodes[elements[i].emap[9]].x +
       4 * w * nodes[elements[i].emap[9]].x - nodes[elements[i].emap[2]].x +
       4 * v * nodes[elements[i].emap[2]].x + nodes[elements[i].emap[3]].x -
       4 * w * nodes[elements[i].emap[3]].x +
       4 * t * nodes[elements[i].emap[5]].x -
       4 * t * nodes[elements[i].emap[6]].x +
       4 * u * nodes[elements[i].emap[7]].x -
       4 * u * nodes[elements[i].emap[8]].x) *
          ((-1 + 4 * t) * nodes[elements[i].emap[0]].z +
           4 * (u * nodes[elements[i].emap[4]].z +
                v * nodes[elements[i].emap[5]].z +
                w * nodes[elements[i].emap[6]].z)) -
      ((-1 + 4 * t) * nodes[elements[i].emap[0]].x -
       4 * v * nodes[elements[i].emap[9]].x + nodes[elements[i].emap[3]].x -
       4 * w * nodes[elements[i].emap[3]].x +
       4 * u * nodes[elements[i].emap[4]].x +
       4 * v * nodes[elements[i].emap[5]].x -
       4 * t * nodes[elements[i].emap[6]].x +
       4 * w * nodes[elements[i].emap[6]].x -
       4 * u * nodes[elements[i].emap[8]].x) *
          (4 * w * nodes[elements[i].emap[9]].z - nodes[elements[i].emap[2]].z +
           4 * v * nodes[elements[i].emap[2]].z +
           4 * t * nodes[elements[i].emap[5]].z +
           4 * u * nodes[elements[i].emap[7]].z) +
      ((-1 + 4 * t) * nodes[elements[i].emap[0]].x -
       4 * w * nodes[elements[i].emap[9]].x + nodes[elements[i].emap[2]].x -
       4 * v * nodes[elements[i].emap[2]].x +
       4 * u * nodes[elements[i].emap[4]].x -
       4 * t * nodes[elements[i].emap[5]].x +
       4 * v * nodes[elements[i].emap[5]].x +
       4 * w * nodes[elements[i].emap[6]].x -
       4 * u * nodes[elements[i].emap[7]].x) *
          (4 * v * nodes[elements[i].emap[9]].z - nodes[elements[i].emap[3]].z +
           4 * w * nodes[elements[i].emap[3]].z +
           4 * t * nodes[elements[i].emap[6]].z +
           4 * u * nodes[elements[i].emap[8]].z);

  jac[1][3] =
      (-4 * w * nodes[elements[i].emap[9]].x +
       4 * v * (nodes[elements[i].emap[9]].x - nodes[elements[i].emap[2]].x) +
       nodes[elements[i].emap[2]].x - nodes[elements[i].emap[3]].x +
       4 * w * nodes[elements[i].emap[3]].x -
       4 * t * nodes[elements[i].emap[5]].x +
       4 * t * nodes[elements[i].emap[6]].x -
       4 * u * nodes[elements[i].emap[7]].x +
       4 * u * nodes[elements[i].emap[8]].x) *
          ((-1 + 4 * t) * nodes[elements[i].emap[0]].y +
           4 * (u * nodes[elements[i].emap[4]].y +
                v * nodes[elements[i].emap[5]].y +
                w * nodes[elements[i].emap[6]].y)) +
      ((-1 + 4 * t) * nodes[elements[i].emap[0]].x -
       4 * v * nodes[elements[i].emap[9]].x + nodes[elements[i].emap[3]].x -
       4 * w * nodes[elements[i].emap[3]].x +
       4 * u * nodes[elements[i].emap[4]].x +
       4 * v * nodes[elements[i].emap[5]].x -
       4 * t * nodes[elements[i].emap[6]].x +
       4 * w * nodes[elements[i].emap[6]].x -
       4 * u * nodes[elements[i].emap[8]].x) *
          (4 * w * nodes[elements[i].emap[9]].y - nodes[elements[i].emap[2]].y +
           4 * v * nodes[elements[i].emap[2]].y +
           4 * t * nodes[elements[i].emap[5]].y +
           4 * u * nodes[elements[i].emap[7]].y) -
      ((-1 + 4 * t) * nodes[elements[i].emap[0]].x -
       4 * w * nodes[elements[i].emap[9]].x + nodes[elements[i].emap[2]].x -
       4 * v * nodes[elements[i].emap[2]].x +
       4 * u * nodes[elements[i].emap[4]].x -
       4 * t * nodes[elements[i].emap[5]].x +
       4 * v * nodes[elements[i].emap[5]].x +
       4 * w * nodes[elements[i].emap[6]].x -
       4 * u * nodes[elements[i].emap[7]].x) *
          (4 * v * nodes[elements[i].emap[9]].y - nodes[elements[i].emap[3]].y +
           4 * w * nodes[elements[i].emap[3]].y +
           4 * t * nodes[elements[i].emap[6]].y +
           4 * u * nodes[elements[i].emap[8]].y);

  jac[2][0] =
      (((-1 + 4 * u) * nodes[elements[i].emap[1]].x +
        4 * (t * nodes[elements[i].emap[4]].x +
             v * nodes[elements[i].emap[7]].x +
             w * nodes[elements[i].emap[8]].x)) *
           (4 * v * nodes[elements[i].emap[9]].y -
            nodes[elements[i].emap[3]].y +
            4 * w * nodes[elements[i].emap[3]].y +
            4 * t * nodes[elements[i].emap[6]].y +
            4 * u * nodes[elements[i].emap[8]].y) -
       (4 * v * nodes[elements[i].emap[9]].x - nodes[elements[i].emap[3]].x +
        4 * w * nodes[elements[i].emap[3]].x +
        4 * t * nodes[elements[i].emap[6]].x +
        4 * u * nodes[elements[i].emap[8]].x) *
           ((-1 + 4 * u) * nodes[elements[i].emap[1]].y +
            4 * (t * nodes[elements[i].emap[4]].y +
                 v * nodes[elements[i].emap[7]].y +
                 w * nodes[elements[i].emap[8]].y))) *
          ((-1 + 4 * t) * nodes[elements[i].emap[0]].z +
           4 * (u * nodes[elements[i].emap[4]].z +
                v * nodes[elements[i].emap[5]].z +
                w * nodes[elements[i].emap[6]].z)) +
      (-(((-1 + 4 * u) * nodes[elements[i].emap[1]].x +
          4 * (t * nodes[elements[i].emap[4]].x +
               v * nodes[elements[i].emap[7]].x +
               w * nodes[elements[i].emap[8]].x)) *
         ((-1 + 4 * t) * nodes[elements[i].emap[0]].y +
          4 * (u * nodes[elements[i].emap[4]].y +
               v * nodes[elements[i].emap[5]].y +
               w * nodes[elements[i].emap[6]].y))) +
       ((-1 + 4 * t) * nodes[elements[i].emap[0]].x +
        4 * (u * nodes[elements[i].emap[4]].x +
             v * nodes[elements[i].emap[5]].x +
             w * nodes[elements[i].emap[6]].x)) *
           ((-1 + 4 * u) * nodes[elements[i].emap[1]].y +
            4 * (t * nodes[elements[i].emap[4]].y +
                 v * nodes[elements[i].emap[7]].y +
                 w * nodes[elements[i].emap[8]].y))) *
          (4 * v * nodes[elements[i].emap[9]].z - nodes[elements[i].emap[3]].z +
           4 * w * nodes[elements[i].emap[3]].z +
           4 * t * nodes[elements[i].emap[6]].z +
           4 * u * nodes[elements[i].emap[8]].z) -
      (-((4 * v * nodes[elements[i].emap[9]].x - nodes[elements[i].emap[3]].x +
          4 * w * nodes[elements[i].emap[3]].x +
          4 * t * nodes[elements[i].emap[6]].x +
          4 * u * nodes[elements[i].emap[8]].x) *
         ((-1 + 4 * t) * nodes[elements[i].emap[0]].y +
          4 * (u * nodes[elements[i].emap[4]].y +
               v * nodes[elements[i].emap[5]].y +
               w * nodes[elements[i].emap[6]].y))) +
       ((-1 + 4 * t) * nodes[elements[i].emap[0]].x +
        4 * (u * nodes[elements[i].emap[4]].x +
             v * nodes[elements[i].emap[5]].x +
             w * nodes[elements[i].emap[6]].x)) *
           (4 * v * nodes[elements[i].emap[9]].y -
            nodes[elements[i].emap[3]].y +
            4 * w * nodes[elements[i].emap[3]].y +
            4 * t * nodes[elements[i].emap[6]].y +
            4 * u * nodes[elements[i].emap[8]].y)) *
          ((-1 + 4 * u) * nodes[elements[i].emap[1]].z +
           4 * (t * nodes[elements[i].emap[4]].z +
                v * nodes[elements[i].emap[7]].z +
                w * nodes[elements[i].emap[8]].z));

  jac[2][1] =
      (-4 * v * nodes[elements[i].emap[9]].y - nodes[elements[i].emap[1]].y +
       4 * u * nodes[elements[i].emap[1]].y + nodes[elements[i].emap[3]].y -
       4 * w * nodes[elements[i].emap[3]].y +
       4 * t * nodes[elements[i].emap[4]].y -
       4 * t * nodes[elements[i].emap[6]].y +
       4 * v * nodes[elements[i].emap[7]].y -
       4 * u * nodes[elements[i].emap[8]].y +
       4 * w * nodes[elements[i].emap[8]].y) *
          ((-1 + 4 * t) * nodes[elements[i].emap[0]].z +
           4 * (u * nodes[elements[i].emap[4]].z +
                v * nodes[elements[i].emap[5]].z +
                w * nodes[elements[i].emap[6]].z)) +
      ((-1 + 4 * t) * nodes[elements[i].emap[0]].y +
       nodes[elements[i].emap[1]].y - 4 * u * nodes[elements[i].emap[1]].y +
       4 * (-(t * nodes[elements[i].emap[4]].y) +
            u * nodes[elements[i].emap[4]].y +
            v * nodes[elements[i].emap[5]].y +
            w * nodes[elements[i].emap[6]].y -
            v * nodes[elements[i].emap[7]].y -
            w * nodes[elements[i].emap[8]].y)) *
          (4 * v * nodes[elements[i].emap[9]].z - nodes[elements[i].emap[3]].z +
           4 * w * nodes[elements[i].emap[3]].z +
           4 * t * nodes[elements[i].emap[6]].z +
           4 * u * nodes[elements[i].emap[8]].z) -
      ((-1 + 4 * t) * nodes[elements[i].emap[0]].y -
       4 * v * nodes[elements[i].emap[9]].y + nodes[elements[i].emap[3]].y -
       4 * w * nodes[elements[i].emap[3]].y +
       4 * u * nodes[elements[i].emap[4]].y +
       4 * v * nodes[elements[i].emap[5]].y -
       4 * t * nodes[elements[i].emap[6]].y +
       4 * w * nodes[elements[i].emap[6]].y -
       4 * u * nodes[elements[i].emap[8]].y) *
          ((-1 + 4 * u) * nodes[elements[i].emap[1]].z +
           4 * (t * nodes[elements[i].emap[4]].z +
                v * nodes[elements[i].emap[7]].z +
                w * nodes[elements[i].emap[8]].z));

  jac[2][2] =
      (nodes[elements[i].emap[1]].x - 4 * u * nodes[elements[i].emap[1]].x -
       nodes[elements[i].emap[3]].x + 4 * w * nodes[elements[i].emap[3]].x -
       4 * t * nodes[elements[i].emap[4]].x +
       4 * t * nodes[elements[i].emap[6]].x +
       4 * v * (nodes[elements[i].emap[9]].x - nodes[elements[i].emap[7]].x) +
       4 * u * nodes[elements[i].emap[8]].x -
       4 * w * nodes[elements[i].emap[8]].x) *
          ((-1 + 4 * t) * nodes[elements[i].emap[0]].z +
           4 * (u * nodes[elements[i].emap[4]].z +
                v * nodes[elements[i].emap[5]].z +
                w * nodes[elements[i].emap[6]].z)) -
      ((-1 + 4 * t) * nodes[elements[i].emap[0]].x +
       nodes[elements[i].emap[1]].x - 4 * u * nodes[elements[i].emap[1]].x +
       4 * (-(t * nodes[elements[i].emap[4]].x) +
            u * nodes[elements[i].emap[4]].x +
            v * nodes[elements[i].emap[5]].x +
            w * nodes[elements[i].emap[6]].x -
            v * nodes[elements[i].emap[7]].x -
            w * nodes[elements[i].emap[8]].x)) *
          (4 * v * nodes[elements[i].emap[9]].z - nodes[elements[i].emap[3]].z +
           4 * w * nodes[elements[i].emap[3]].z +
           4 * t * nodes[elements[i].emap[6]].z +
           4 * u * nodes[elements[i].emap[8]].z) +
      ((-1 + 4 * t) * nodes[elements[i].emap[0]].x -
       4 * v * nodes[elements[i].emap[9]].x + nodes[elements[i].emap[3]].x -
       4 * w * nodes[elements[i].emap[3]].x +
       4 * u * nodes[elements[i].emap[4]].x +
       4 * v * nodes[elements[i].emap[5]].x -
       4 * t * nodes[elements[i].emap[6]].x +
       4 * w * nodes[elements[i].emap[6]].x -
       4 * u * nodes[elements[i].emap[8]].x) *
          ((-1 + 4 * u) * nodes[elements[i].emap[1]].z +
           4 * (t * nodes[elements[i].emap[4]].z +
                v * nodes[elements[i].emap[7]].z +
                w * nodes[elements[i].emap[8]].z));

  jac[2][3] =
      (-4 * v * nodes[elements[i].emap[9]].x - nodes[elements[i].emap[1]].x +
       4 * u * nodes[elements[i].emap[1]].x + nodes[elements[i].emap[3]].x -
       4 * w * nodes[elements[i].emap[3]].x +
       4 * t * nodes[elements[i].emap[4]].x -
       4 * t * nodes[elements[i].emap[6]].x +
       4 * v * nodes[elements[i].emap[7]].x -
       4 * u * nodes[elements[i].emap[8]].x +
       4 * w * nodes[elements[i].emap[8]].x) *
          ((-1 + 4 * t) * nodes[elements[i].emap[0]].y +
           4 * (u * nodes[elements[i].emap[4]].y +
                v * nodes[elements[i].emap[5]].y +
                w * nodes[elements[i].emap[6]].y)) +
      ((-1 + 4 * t) * nodes[elements[i].emap[0]].x +
       nodes[elements[i].emap[1]].x - 4 * u * nodes[elements[i].emap[1]].x +
       4 * (-(t * nodes[elements[i].emap[4]].x) +
            u * nodes[elements[i].emap[4]].x +
            v * nodes[elements[i].emap[5]].x +
            w * nodes[elements[i].emap[6]].x -
            v * nodes[elements[i].emap[7]].x -
            w * nodes[elements[i].emap[8]].x)) *
          (4 * v * nodes[elements[i].emap[9]].y - nodes[elements[i].emap[3]].y +
           4 * w * nodes[elements[i].emap[3]].y +
           4 * t * nodes[elements[i].emap[6]].y +
           4 * u * nodes[elements[i].emap[8]].y) -
      ((-1 + 4 * t) * nodes[elements[i].emap[0]].x -
       4 * v * nodes[elements[i].emap[9]].x + nodes[elements[i].emap[3]].x -
       4 * w * nodes[elements[i].emap[3]].x +
       4 * u * nodes[elements[i].emap[4]].x +
       4 * v * nodes[elements[i].emap[5]].x -
       4 * t * nodes[elements[i].emap[6]].x +
       4 * w * nodes[elements[i].emap[6]].x -
       4 * u * nodes[elements[i].emap[8]].x) *
          ((-1 + 4 * u) * nodes[elements[i].emap[1]].y +
           4 * (t * nodes[elements[i].emap[4]].y +
                v * nodes[elements[i].emap[7]].y +
                w * nodes[elements[i].emap[8]].y));

  jac[3][0] =
      -((((-1 + 4 * u) * nodes[elements[i].emap[1]].x +
          4 * (t * nodes[elements[i].emap[4]].x +
               v * nodes[elements[i].emap[7]].x +
               w * nodes[elements[i].emap[8]].x)) *
             (4 * w * nodes[elements[i].emap[9]].y -
              nodes[elements[i].emap[2]].y +
              4 * v * nodes[elements[i].emap[2]].y +
              4 * t * nodes[elements[i].emap[5]].y +
              4 * u * nodes[elements[i].emap[7]].y) -
         (4 * w * nodes[elements[i].emap[9]].x - nodes[elements[i].emap[2]].x +
          4 * v * nodes[elements[i].emap[2]].x +
          4 * t * nodes[elements[i].emap[5]].x +
          4 * u * nodes[elements[i].emap[7]].x) *
             ((-1 + 4 * u) * nodes[elements[i].emap[1]].y +
              4 * (t * nodes[elements[i].emap[4]].y +
                   v * nodes[elements[i].emap[7]].y +
                   w * nodes[elements[i].emap[8]].y))) *
        ((-1 + 4 * t) * nodes[elements[i].emap[0]].z +
         4 * (u * nodes[elements[i].emap[4]].z +
              v * nodes[elements[i].emap[5]].z +
              w * nodes[elements[i].emap[6]].z))) -
      (-(((-1 + 4 * u) * nodes[elements[i].emap[1]].x +
          4 * (t * nodes[elements[i].emap[4]].x +
               v * nodes[elements[i].emap[7]].x +
               w * nodes[elements[i].emap[8]].x)) *
         ((-1 + 4 * t) * nodes[elements[i].emap[0]].y +
          4 * (u * nodes[elements[i].emap[4]].y +
               v * nodes[elements[i].emap[5]].y +
               w * nodes[elements[i].emap[6]].y))) +
       ((-1 + 4 * t) * nodes[elements[i].emap[0]].x +
        4 * (u * nodes[elements[i].emap[4]].x +
             v * nodes[elements[i].emap[5]].x +
             w * nodes[elements[i].emap[6]].x)) *
           ((-1 + 4 * u) * nodes[elements[i].emap[1]].y +
            4 * (t * nodes[elements[i].emap[4]].y +
                 v * nodes[elements[i].emap[7]].y +
                 w * nodes[elements[i].emap[8]].y))) *
          (4 * w * nodes[elements[i].emap[9]].z - nodes[elements[i].emap[2]].z +
           4 * v * nodes[elements[i].emap[2]].z +
           4 * t * nodes[elements[i].emap[5]].z +
           4 * u * nodes[elements[i].emap[7]].z) +
      (-((4 * w * nodes[elements[i].emap[9]].x - nodes[elements[i].emap[2]].x +
          4 * v * nodes[elements[i].emap[2]].x +
          4 * t * nodes[elements[i].emap[5]].x +
          4 * u * nodes[elements[i].emap[7]].x) *
         ((-1 + 4 * t) * nodes[elements[i].emap[0]].y +
          4 * (u * nodes[elements[i].emap[4]].y +
               v * nodes[elements[i].emap[5]].y +
               w * nodes[elements[i].emap[6]].y))) +
       ((-1 + 4 * t) * nodes[elements[i].emap[0]].x +
        4 * (u * nodes[elements[i].emap[4]].x +
             v * nodes[elements[i].emap[5]].x +
             w * nodes[elements[i].emap[6]].x)) *
           (4 * w * nodes[elements[i].emap[9]].y -
            nodes[elements[i].emap[2]].y +
            4 * v * nodes[elements[i].emap[2]].y +
            4 * t * nodes[elements[i].emap[5]].y +
            4 * u * nodes[elements[i].emap[7]].y)) *
          ((-1 + 4 * u) * nodes[elements[i].emap[1]].z +
           4 * (t * nodes[elements[i].emap[4]].z +
                v * nodes[elements[i].emap[7]].z +
                w * nodes[elements[i].emap[8]].z));

  jac[3][1] =
      (nodes[elements[i].emap[1]].y - 4 * u * nodes[elements[i].emap[1]].y -
       nodes[elements[i].emap[2]].y + 4 * v * nodes[elements[i].emap[2]].y -
       4 * t * nodes[elements[i].emap[4]].y +
       4 * t * nodes[elements[i].emap[5]].y +
       4 * u * nodes[elements[i].emap[7]].y -
       4 * v * nodes[elements[i].emap[7]].y +
       4 * w * (nodes[elements[i].emap[9]].y - nodes[elements[i].emap[8]].y)) *
          ((-1 + 4 * t) * nodes[elements[i].emap[0]].z +
           4 * (u * nodes[elements[i].emap[4]].z +
                v * nodes[elements[i].emap[5]].z +
                w * nodes[elements[i].emap[6]].z)) -
      ((-1 + 4 * t) * nodes[elements[i].emap[0]].y +
       nodes[elements[i].emap[1]].y - 4 * u * nodes[elements[i].emap[1]].y +
       4 * (-(t * nodes[elements[i].emap[4]].y) +
            u * nodes[elements[i].emap[4]].y +
            v * nodes[elements[i].emap[5]].y +
            w * nodes[elements[i].emap[6]].y -
            v * nodes[elements[i].emap[7]].y -
            w * nodes[elements[i].emap[8]].y)) *
          (4 * w * nodes[elements[i].emap[9]].z - nodes[elements[i].emap[2]].z +
           4 * v * nodes[elements[i].emap[2]].z +
           4 * t * nodes[elements[i].emap[5]].z +
           4 * u * nodes[elements[i].emap[7]].z) +
      ((-1 + 4 * t) * nodes[elements[i].emap[0]].y -
       4 * w * nodes[elements[i].emap[9]].y + nodes[elements[i].emap[2]].y -
       4 * v * nodes[elements[i].emap[2]].y +
       4 * u * nodes[elements[i].emap[4]].y -
       4 * t * nodes[elements[i].emap[5]].y +
       4 * v * nodes[elements[i].emap[5]].y +
       4 * w * nodes[elements[i].emap[6]].y -
       4 * u * nodes[elements[i].emap[7]].y) *
          ((-1 + 4 * u) * nodes[elements[i].emap[1]].z +
           4 * (t * nodes[elements[i].emap[4]].z +
                v * nodes[elements[i].emap[7]].z +
                w * nodes[elements[i].emap[8]].z));

  jac[3][2] =
      (-4 * w * nodes[elements[i].emap[9]].x - nodes[elements[i].emap[1]].x +
       4 * u * nodes[elements[i].emap[1]].x + nodes[elements[i].emap[2]].x -
       4 * v * nodes[elements[i].emap[2]].x +
       4 * t * nodes[elements[i].emap[4]].x -
       4 * t * nodes[elements[i].emap[5]].x -
       4 * u * nodes[elements[i].emap[7]].x +
       4 * v * nodes[elements[i].emap[7]].x +
       4 * w * nodes[elements[i].emap[8]].x) *
          ((-1 + 4 * t) * nodes[elements[i].emap[0]].z +
           4 * (u * nodes[elements[i].emap[4]].z +
                v * nodes[elements[i].emap[5]].z +
                w * nodes[elements[i].emap[6]].z)) +
      ((-1 + 4 * t) * nodes[elements[i].emap[0]].x +
       nodes[elements[i].emap[1]].x - 4 * u * nodes[elements[i].emap[1]].x +
       4 * (-(t * nodes[elements[i].emap[4]].x) +
            u * nodes[elements[i].emap[4]].x +
            v * nodes[elements[i].emap[5]].x +
            w * nodes[elements[i].emap[6]].x -
            v * nodes[elements[i].emap[7]].x -
            w * nodes[elements[i].emap[8]].x)) *
          (4 * w * nodes[elements[i].emap[9]].z - nodes[elements[i].emap[2]].z +
           4 * v * nodes[elements[i].emap[2]].z +
           4 * t * nodes[elements[i].emap[5]].z +
           4 * u * nodes[elements[i].emap[7]].z) -
      ((-1 + 4 * t) * nodes[elements[i].emap[0]].x -
       4 * w * nodes[elements[i].emap[9]].x + nodes[elements[i].emap[2]].x -
       4 * v * nodes[elements[i].emap[2]].x +
       4 * u * nodes[elements[i].emap[4]].x -
       4 * t * nodes[elements[i].emap[5]].x +
       4 * v * nodes[elements[i].emap[5]].x +
       4 * w * nodes[elements[i].emap[6]].x -
       4 * u * nodes[elements[i].emap[7]].x) *
          ((-1 + 4 * u) * nodes[elements[i].emap[1]].z +
           4 * (t * nodes[elements[i].emap[4]].z +
                v * nodes[elements[i].emap[7]].z +
                w * nodes[elements[i].emap[8]].z));

  jac[3][3] =
      (nodes[elements[i].emap[1]].x - 4 * u * nodes[elements[i].emap[1]].x -
       nodes[elements[i].emap[2]].x + 4 * v * nodes[elements[i].emap[2]].x -
       4 * t * nodes[elements[i].emap[4]].x +
       4 * t * nodes[elements[i].emap[5]].x +
       4 * u * nodes[elements[i].emap[7]].x -
       4 * v * nodes[elements[i].emap[7]].x +
       4 * w * (nodes[elements[i].emap[9]].x - nodes[elements[i].emap[8]].x)) *
          ((-1 + 4 * t) * nodes[elements[i].emap[0]].y +
           4 * (u * nodes[elements[i].emap[4]].y +
                v * nodes[elements[i].emap[5]].y +
                w * nodes[elements[i].emap[6]].y)) -
      ((-1 + 4 * t) * nodes[elements[i].emap[0]].x +
       nodes[elements[i].emap[1]].x - 4 * u * nodes[elements[i].emap[1]].x +
       4 * (-(t * nodes[elements[i].emap[4]].x) +
            u * nodes[elements[i].emap[4]].x +
            v * nodes[elements[i].emap[5]].x +
            w * nodes[elements[i].emap[6]].x -
            v * nodes[elements[i].emap[7]].x -
            w * nodes[elements[i].emap[8]].x)) *
          (4 * w * nodes[elements[i].emap[9]].y - nodes[elements[i].emap[2]].y +
           4 * v * nodes[elements[i].emap[2]].y +
           4 * t * nodes[elements[i].emap[5]].y +
           4 * u * nodes[elements[i].emap[7]].y) +
      ((-1 + 4 * t) * nodes[elements[i].emap[0]].x -
       4 * w * nodes[elements[i].emap[9]].x + nodes[elements[i].emap[2]].x -
       4 * v * nodes[elements[i].emap[2]].x +
       4 * u * nodes[elements[i].emap[4]].x -
       4 * t * nodes[elements[i].emap[5]].x +
       4 * v * nodes[elements[i].emap[5]].x +
       4 * w * nodes[elements[i].emap[6]].x -
       4 * u * nodes[elements[i].emap[7]].x) *
          ((-1 + 4 * u) * nodes[elements[i].emap[1]].y +
           4 * (t * nodes[elements[i].emap[4]].y +
                v * nodes[elements[i].emap[7]].y +
                w * nodes[elements[i].emap[8]].y));
}
