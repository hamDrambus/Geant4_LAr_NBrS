#include <field_drift/FieldElmerMap.hh>

FieldElmerMap::FieldElmerMap() :
    nElements(-1),
    lastElement(-1),
    nNodes(0),
    nMaterials(0),
    hasBoundingBox(false),
    checkMultipleElement(false),
    warning(false),
    m_octree(nullptr),
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
  Node newNode;
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
  Element newElement;
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
      Material mat;
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
  InitializeTetrahedralTree();

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
    Element& elem = elements[i];
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

bool FieldElmerMap::InitializeTetrahedralTree(void)
{
  // Do not proceed if not properly initialised.
  if (!ready) {
    std::cerr << m_className << "::InitializeTetrahedralTree:\n";
    std::cerr << "    Field map not yet initialised.\n";
    std::cerr << "    Bounding boxes of elements cannot be computed.\n";
    return false;
  }

  if (debug) {
    std::cout << m_className << "::InitializeTetrahedralTree:\n"
        << "    About to initialize the tetrahedral tree.\n";
  }

  if (nodes.empty()) {
    std::cerr << m_className << "::InitializeTetrahedralTree: Empty mesh.\n";
    return false;
  }

  // Determine the bounding box
  double xmin = nodes.front().x;
  double ymin = nodes.front().y;
  double zmin = nodes.front().z;
  double xmax = xmin;
  double ymax = ymin;
  double zmax = zmin;
  for (const auto& node : nodes) {
    xmin = std::min(xmin, node.x);
    xmax = std::max(xmax, node.x);
    ymin = std::min(ymin, node.y);
    ymax = std::max(ymax, node.y);
    zmin = std::min(zmin, node.z);
    zmax = std::max(zmax, node.z);
  }

  if (debug) {
    std::cout << "    Bounding box:\n"
        << std::scientific << "\tx: " << xmin << " -> " << xmax << "\n"
        << std::scientific << "\ty: " << ymin << " -> " << ymax << "\n"
        << std::scientific << "\tz: " << zmin << " -> " << zmax << "\n";
  }

  const double hx = 0.5 * (xmax - xmin);
  const double hy = 0.5 * (ymax - ymin);
  const double hz = 0.5 * (zmax - zmin);
  m_octree.reset(new TetrahedralTree(G4ThreeVector(xmin + hx, ymin + hy, zmin + hz),
      G4ThreeVector(hx, hy, hz)));

  if (debug)
    std::cout << "    Tree instantiated.\n";

  // Insert all mesh nodes in the tree
  for (unsigned int i = 0; i < nodes.size(); i++) {
    const Node& n = nodes[i];
    m_octree->InsertMeshNode(G4ThreeVector(n.x, n.y, n.z), i);
  }

  if (debug)
    std::cout << "    Tree nodes initialized successfully.\n";

  // Insert all mesh elements (tetrahedrons) in the tree
  for (unsigned int i = 0; i < elements.size(); ++i) {
    const Element& e = elements[i];
    const double bb[6] = {e.xmin, e.ymin, e.zmin,
                          e.xmax, e.ymax, e.zmax};
    m_octree->InsertMeshElement(bb, i);
  }

  m_octree->Finalize();

  std::cout << m_className << "::InitializeTetrahedralTree: Success.\n";
  return true;
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
  const Element& element = elements[imap];
  std::array<double, 10> v;
  for (size_t i = 0; i < 10; ++i)
    v[i] = nodes[element.emap[i]].v;
  volt = Potential13(v, {t1, t2, t3, t4});
  Field13(v, {t1, t2, t3, t4}, jac, det, ex, ey, ez);

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

double FieldElmerMap::Potential13(const std::array<double, 10>& v, const std::array<double, 4>& t)
{
  double sum = 0.;
  for (size_t i = 0; i < 4; ++i) {
    sum += v[i] * t[i] * (t[i] - 0.5);
  }
  sum *= 2;
  sum += 4 * (v[4] * t[0] * t[1] + v[5] * t[0] * t[2] + v[6] * t[0] * t[3] +
              v[7] * t[1] * t[2] + v[8] * t[1] * t[3] + v[9] * t[2] * t[3]);
  return sum;
}

void FieldElmerMap::Field13(const std::array<double, 10>& v, const std::array<double, 4>& t,
                 double jac[4][4], const double det, double& ex, double& ey, double& ez)
{
  std::array<double, 4> g;
  g[0] = v[0] * (t[0] - 0.25) + v[4] * t[1] + v[5] * t[2] + v[6] * t[3];
  g[1] = v[1] * (t[1] - 0.25) + v[4] * t[0] + v[7] * t[2] + v[8] * t[3];
  g[2] = v[2] * (t[2] - 0.25) + v[5] * t[0] + v[7] * t[1] + v[9] * t[3];
  g[3] = v[3] * (t[3] - 0.25) + v[6] * t[0] + v[8] * t[1] + v[9] * t[2];
  std::array<double, 3> f = {0., 0., 0.};
  for (size_t j = 0; j < 4; ++j) {
    for (size_t i = 0; i < 3; ++i) {
      f[i] += g[j] * jac[j][i + 1];
    }
  }
  const double invdet = -4. / det;
  ex = f[0] * invdet;
  ey = f[1] * invdet;
  ez = f[2] * invdet;
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

  // The loop below can be written without lambda function but it will be slower by ~35%.
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

  if (!m_useTetrahedralTree || !m_octree) {
    for (int i = 0; i != nElements; ++i) {
      Element& e = elements[i];
      if (x < e.xmin || x > e.xmax ||
          y < e.ymin || y > e.ymax ||
          z < e.zmin || z > e.zmax)
        continue;
      int res = process_found_element(i);
      if (res >= 0)
        return res;
    }
  } else {
    std::vector<int> tetList = m_octree->GetElementsInBlock(G4ThreeVector(x, y, z));
    for (int j = 0, j_end_ = tetList.size(); j != j_end_; ++j) {
      int i = tetList[j];
      Element& e = elements[i];
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

  // Make a first order approximation.
  Coordinates12(x, y, z, t1, t2, t3, t4, imap);

  // Set tolerance parameter.
  constexpr double f = 0.5;
  if (t1 < -f || t2 < -f || t3 < -f || t4 < -f || t1 > 1 + f || t2 > 1 + f ||
      t3 > 1 + f || t4 > 1 + f) {
    if (debug) {
      std::cout << m_className << "::Coordinates13:\n"
                << "    Linear isoparametric coordinates more than\n"
                << "    f (" << f << ") out of range.\n";
    }
    return 0;
  }

  // (Mine)
  //t1 = std::min(std::max(t1, 0.0), 1.0);
  //t2 = std::min(std::max(t2, 0.0), 1.0);
  //t3 = std::min(std::max(t3, 0.0), 1.0);
  //t4 = std::min(std::max(t4, 0.0), 1.0);
  // Start iteration.
  double td1 = t1, td2 = t2, td3 = t3, td4 = t4;
  const Element & element = elements[imap];
  const Node& n0 = nodes[element.emap[0]];
  const Node& n1 = nodes[element.emap[1]];
  const Node& n2 = nodes[element.emap[2]];
  const Node& n3 = nodes[element.emap[3]];
  const Node& n4 = nodes[element.emap[4]];
  const Node& n5 = nodes[element.emap[5]];
  const Node& n6 = nodes[element.emap[6]];
  const Node& n7 = nodes[element.emap[7]];
  const Node& n8 = nodes[element.emap[8]];
  const Node& n9 = nodes[element.emap[9]];

  // Loop
  bool converged = false;
  double diff[4], corr[4];
  for (int iter = 0; iter < 10; iter++) {
    if (debug) {
      std::cout << m_className << "::Coordinates13:\n";
      std::printf("    Iteration %4u: t = (%15.8f, %15.8f %15.8f %15.8f)\n",
                  iter, td1, td2, td3, td4);
    }
    const double f0 = td1 * (2 * td1 - 1);
    const double f1 = td2 * (2 * td2 - 1);
    const double f2 = td3 * (2 * td3 - 1);
    const double f3 = td4 * (2 * td4 - 1);
    const double f4 = 4 * td1 * td2;
    const double f5 = 4 * td1 * td3;
    const double f6 = 4 * td1 * td4;
    const double f7 = 4 * td2 * td3;
    const double f8 = 4 * td2 * td4;
    const double f9 = 4 * td3 * td4;
    // Re-compute the (x,y,z) position for this coordinate.
    const double xr = n0.x * f0 + n1.x * f1 + n2.x * f2 + n3.x * f3 +
                      n4.x * f4 + n5.x * f5 + n6.x * f6 + n7.x * f7 +
                      n8.x * f8 + n9.x * f9;
    const double yr = n0.y * f0 + n1.y * f1 + n2.y * f2 + n3.y * f3 +
                      n4.y * f4 + n5.y * f5 + n6.y * f6 + n7.y * f7 +
                      n8.y * f8 + n9.y * f9;
    const double zr = n0.z * f0 + n1.z * f1 + n2.z * f2 + n3.z * f3 +
                      n4.z * f4 + n5.z * f5 + n6.z * f6 + n7.z * f7 +
                      n8.z * f8 + n9.z * f9;
    const double sr = td1 + td2 + td3 + td4;

    // Compute the Jacobian.
    Jacobian13(imap, td1, td2, td3, td4, det, jac);
    // Compute the difference vector.
    diff[0] = 1 - sr;
    diff[1] = x - xr;
    diff[2] = y - yr;
    diff[3] = z - zr;
    // Update the estimate.
    const double invdet = 1. / det;
    for (int l = 0; l < 4; ++l) {
      corr[l] = 0;
      for (int k = 0; k < 4; ++k) {
        corr[l] += jac[l][k] * diff[k];
      }
      corr[l] *= invdet;
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
    constexpr double tol = 1.e-5;
    if (fabs(corr[0]) < tol && fabs(corr[1]) < tol && fabs(corr[2]) < tol &&
        fabs(corr[3]) < tol) {
      if (debug) {
        std::cout << m_className << "::Coordinates13: Convergence reached.\n";
      }
      converged = true;
      break;
    }
  }

  // No convergence reached.
  if (!converged) {
    const double xmin = std::min({n0.x, n1.x, n2.x, n3.x});
    const double xmax = std::max({n0.x, n1.x, n2.x, n3.x});
    const double ymin = std::min({n0.y, n1.y, n2.y, n3.y});
    const double ymax = std::max({n0.y, n1.y, n2.y, n3.y});
    const double zmin = std::min({n0.z, n1.z, n2.z, n3.z});
    const double zmax = std::max({n0.z, n1.z, n2.z, n3.z});
    if (x >= xmin && x <= xmax && y >= ymin && y <= ymax && z >= zmin &&
        z <= zmax) {
      std::cout << m_className << "::Coordinates13:\n"
                << "    No convergence achieved "
                << "when refining internal isoparametric coordinates\n"
                << "    at position (" << x << ", " << y << ", " << z
                << ").\n";
      t1 = t2 = t3 = t4 = -1;
      return 1;
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
    // Re-compute the (x,y,z) position for this coordinate.
    const double f0 = td1 * (2 * td1 - 1);
    const double f1 = td2 * (2 * td2 - 1);
    const double f2 = td3 * (2 * td3 - 1);
    const double f3 = td4 * (2 * td4 - 1);
    const double f4 = 4 * td1 * td2;
    const double f5 = 4 * td1 * td3;
    const double f6 = 4 * td1 * td4;
    const double f7 = 4 * td2 * td3;
    const double f8 = 4 * td2 * td4;
    const double f9 = 4 * td3 * td4;
    double xr = n0.x * f0 + n1.x * f1 + n2.x * f2 + n3.x * f3 + n4.x * f4 +
                n5.x * f5 + n6.x * f6 + n7.x * f7 + n8.x * f8 + n9.x * f9;
    double yr = n0.y * f0 + n1.y * f1 + n2.y * f2 + n3.y * f3 + n4.y * f4 +
                n5.y * f5 + n6.y * f6 + n7.y * f7 + n8.y * f8 + n9.y * f9;
    double zr = n0.z * f0 + n1.z * f1 + n2.z * f2 + n3.z * f3 + n4.z * f4 +
                n5.z * f5 + n6.z * f6 + n7.z * f7 + n8.z * f8 + n9.z * f9;
    double sr = td1 + td2 + td3 + td4;
    std::cout << "    Position requested:     (" << x << ", " << y << ", " << z
              << ")\n";
    std::cout << "    Reconstructed:          (" << xr << ", " << yr << ", "
              << zr << ")\n";
    std::cout << "    Difference:             (" << x - xr << ", " << y - yr
              << ", " << z - zr << ")\n";
    std::cout << "    Checksum - 1:           " << sr - 1 << "\n";
  }

  // Success
  return 0;
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
  const Element& element = elements[imap];
  const Node& n0 = nodes[element.emap[0]];
  const Node& n1 = nodes[element.emap[1]];
  const Node& n2 = nodes[element.emap[2]];
  const Node& n3 = nodes[element.emap[3]];
  // Compute tetrahedral coordinates.
  const double f1x =
      (n2.y - n1.y) * (n3.z - n1.z) - (n3.y - n1.y) * (n2.z - n1.z);
  const double f1y =
      (n2.z - n1.z) * (n3.x - n1.x) - (n3.z - n1.z) * (n2.x - n1.x);
  const double f1z =
      (n2.x - n1.x) * (n3.y - n1.y) - (n3.x - n1.x) * (n2.y - n1.y);
  t1 = (x - n1.x) * f1x + (y - n1.y) * f1y + (z - n1.z) * f1z;
  t1 = t1 / ((n0.x - n1.x) * f1x + (n0.y - n1.y) * f1y + (n0.z - n1.z) * f1z);
  const double f2x =
      (n0.y - n2.y) * (n3.z - n2.z) - (n3.y - n2.y) * (n0.z - n2.z);
  const double f2y =
      (n0.z - n2.z) * (n3.x - n2.x) - (n3.z - n2.z) * (n0.x - n2.x);
  const double f2z =
      (n0.x - n2.x) * (n3.y - n2.y) - (n3.x - n2.x) * (n0.y - n2.y);
  t2 = (x - n2.x) * f2x + (y - n2.y) * f2y + (z - n2.z) * f2z;
  t2 = t2 / ((n1.x - n2.x) * f2x + (n1.y - n2.y) * f2y + (n1.z - n2.z) * f2z);
  const double f3x =
      (n0.y - n3.y) * (n1.z - n3.z) - (n1.y - n3.y) * (n0.z - n3.z);
  const double f3y =
      (n0.z - n3.z) * (n1.x - n3.x) - (n1.z - n3.z) * (n0.x - n3.x);
  const double f3z =
      (n0.x - n3.x) * (n1.y - n3.y) - (n1.x - n3.x) * (n0.y - n3.y);
  t3 = (x - n3.x) * f3x + (y - n3.y) * f3y + (z - n3.z) * f3z;
  t3 = t3 / ((n2.x - n3.x) * f3x + (n2.y - n3.y) * f3y + (n2.z - n3.z) * f3z);
  const double f4x =
      (n2.y - n0.y) * (n1.z - n0.z) - (n1.y - n0.y) * (n2.z - n0.z);
  const double f4y =
      (n2.z - n0.z) * (n1.x - n0.x) - (n1.z - n0.z) * (n2.x - n0.x);
  const double f4z =
      (n2.x - n0.x) * (n1.y - n0.y) - (n1.x - n0.x) * (n2.y - n0.y);
  t4 = (x - n0.x) * f4x + (y - n0.y) * f4y + (z - n0.z) * f4z;
  t4 = t4 / ((n3.x - n0.x) * f4x + (n3.y - n0.y) * f4y + (n3.z - n0.z) * f4z);


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

  const Element& element = elements[i];

  const Node& n0 = nodes[element.emap[0]];
  const Node& n1 = nodes[element.emap[1]];
  const Node& n2 = nodes[element.emap[2]];
  const Node& n3 = nodes[element.emap[3]];
  const Node& n4 = nodes[element.emap[4]];
  const Node& n5 = nodes[element.emap[5]];
  const Node& n6 = nodes[element.emap[6]];
  const Node& n7 = nodes[element.emap[7]];
  const Node& n8 = nodes[element.emap[8]];
  const Node& n9 = nodes[element.emap[9]];

  const double tx = 4 * ((-0.25 + t) * n0.x + u * n4.x + v * n5.x + w * n6.x);
  const double ty = 4 * ((-0.25 + t) * n0.y + u * n4.y + v * n5.y + w * n6.y);
  const double tz = 4 * ((-0.25 + t) * n0.z + u * n4.z + v * n5.z + w * n6.z);

  const double ux = 4 * ((-0.25 + u) * n1.x + t * n4.x + v * n7.x + w * n8.x);
  const double uy = 4 * ((-0.25 + u) * n1.y + t * n4.y + v * n7.y + w * n8.y);
  const double uz = 4 * ((-0.25 + u) * n1.z + t * n4.z + v * n7.z + w * n8.z);

  const double vx = 4 * ((-0.25 + v) * n2.x + t * n5.x + u * n7.x + w * n9.x);
  const double vy = 4 * ((-0.25 + v) * n2.y + t * n5.y + u * n7.y + w * n9.y);
  const double vz = 4 * ((-0.25 + v) * n2.z + t * n5.z + u * n7.z + w * n9.z);

  const double wx = 4 * ((-0.25 + w) * n3.x + t * n6.x + u * n8.x + v * n9.x);
  const double wy = 4 * ((-0.25 + w) * n3.y + t * n6.y + u * n8.y + v * n9.y);
  const double wz = 4 * ((-0.25 + w) * n3.z + t * n6.z + u * n8.z + v * n9.z);

  const double ax = ux - wx;
  const double ay = uy - wy;

  const double bx = ux - vx;
  const double by = uy - vy;

  const double cx = vx - wx;
  const double cy = vy - wy;

  const double dx = tx - wx;
  const double dy = ty - wy;

  const double ex = tx - vx;
  const double ey = ty - vy;

  const double fx = tx - ux;
  const double fy = ty - uy;

  // Determinant of the quadrilateral serendipity Jacobian
  det = (-ax * vy + bx * wy + cx * uy) * tz -
        (-ax * ty - fx * wy + dx * uy) * vz +
        (-bx * ty - fx * vy + ex * uy) * wz +
        (-cx * ty + dx * vy - ex * wy) * uz;

  const double tu = tx * uy - ux * ty;
  const double tv = tx * vy - vx * ty;
  const double tw = tx * wy - wx * ty;
  const double uv = ux * vy - vx * uy;
  const double uw = ux * wy - wx * uy;
  const double vw = vx * wy - wx * vy;

  jac[0][0] = -uw * vz + uv * wz + vw * uz;
  jac[1][0] = -vw * tz + tw * vz - tv * wz;
  jac[2][0] =  uw * tz + tu * wz - tw * uz;
  jac[3][0] = -uv * tz - tu * vz + tv * uz;

  jac[0][1] = -ay * vz + by * wz + cy * uz;
  jac[0][2] =  ax * vz - bx * wz - cx * uz;
  jac[0][3] = -ax * vy + bx * wy + cx * uy;

  jac[1][1] = -cy * tz + dy * vz - ey * wz;
  jac[1][2] =  cx * tz - dx * vz + ex * wz;
  jac[1][3] = -cx * ty + dx * vy - ex * wy;

  jac[2][1] =  ay * tz + fy * wz - dy * uz;
  jac[2][2] = -ax * tz - fx * wz + dx * uz;
  jac[2][3] =  ax * ty + fx * wy - dx * uy;

  jac[3][1] = -by * tz - fy * vz + ey * uz;
  jac[3][2] =  bx * tz + fx * vz - ex * uz;
  jac[3][3] = -bx * ty - fx * vy + ex * uy;
}
