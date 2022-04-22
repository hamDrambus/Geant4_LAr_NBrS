/*  This is copy of Garfield++ ComponentElmer and its base classes
 *  with unnecessary code (for more general field maps) removed.
 *  CERN ROOT (for matrix operations) dependency was also replaced for Boost.
 *  Hence the requirements for Elmer output file formats are the same as
 *  for Garfield++.
 */

#ifndef FIELD_ELMER_MAP_H_
#define FIELD_ELMER_MAP_H_

#include <string>
#include <vector>
#include <deque>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <G4Cache.hh>

#include "DriftMedium.hh"


//TODO: For speeding up element search use new Garfield++'s TetrahedralTree class.
//Check it against my approach. It will possibly need to be optimized as well (remove recursion)
class FieldElmerMap {
public:
  FieldElmerMap();
  FieldElmerMap(std::string header, std::string elist, std::string nlist,
                std::string volt, std::string unit);
  ~FieldElmerMap();
  bool Initialise(std::string header = "mesh.header",
                   std::string elist = "mesh.elements",
                   std::string nlist = "mesh.nodes",
                   std::string volt = "out.result", std::string unit = "mm");

  // Thread-safe event though it is not const
  void ElectricField(const double x, const double y, const double z, double& ex,
                    double& ey, double& ez, DriftMedium* &medium, int& status);
  void ElectricField(const double x, const double y, const double z, double& ex,
                    double& ey, double& ez, double& v, DriftMedium* &medium, int& status);

  // Ranges
  // Calculates x, y, z, V and angular ranges
  void SetRange();
  // Shows x, y, z, V and angular ranges
  void PrintRange() const;

  bool GetBoundingBox(double& xmin, double& ymin, double& zmin,
      double& xmax, double& ymax, double& zmax) const
  {
    if (!ready) return false;
    xmin = mapxmin;
    xmax = mapxmax;
    ymin = mapymin;
    ymax = mapymax;
    zmin = mapzmin;
    zmax = mapzmax;
    return true;
  }

  bool GetVoltageRange(double& vmin, double& vmax) const {
    vmin = mapvmin;
    vmax = mapvmax;
    return true;
  }

  bool IsInBoundingBox(const G4ThreeVector &pos) const {
    return pos.x() >= mapxmin && pos.x() <= mapxmax &&
        pos.y() >= mapymin && pos.y() <= mapymax &&
        pos.z() >= mapzmin && pos.z() <= mapzmax;
  }

  int GetNumberOfElements() const { return nElements; }
  int GetNumberOfMaterials() const { return nMaterials; }
  // Associate a material with a Medium class
  void SetMedium(const int imat, DriftMedium* medium);
  // Returns the medium for a material
  DriftMedium* GetMedium(const unsigned int& i) const;

  // Ready for use?
  bool IsReady() const { return ready; }

  // Options
  void EnableCheckMapIndices()
  {
   checkMultipleElement = true;
   lastElement = -1;
  }
  void DisableCheckMapIndices() { checkMultipleElement = false; }
  void SetRelativeTolerance(double tolerance) { fTolerance = std::fabs(tolerance); }
  double GetRelativeTolerance(void) const { return fTolerance; }
  void EnableDebugging() { debug = true; }
  void DisableDebugging() { debug = false; }

protected:
  std::string m_className;
  bool ready;
  double fTolerance;

  // Grouping of elements into regions for faster search of element corresponding to (x,y,z)
  int nRegions;
  struct region {
    std::vector<std::size_t> indices; // indices of elements inside this region
    double xmin, ymin, zmin, xmax, ymax, zmax; // Bounding box of the region and all contained elements
  };
  std::deque<region> regions;

  // Elements
  int nElements;
  struct element {
   int emap[10]; // Nodes
   int matmap; // Material
   bool degenerate;
   // Bounding box of the element
   double xmin, ymin, zmin, xmax, ymax, zmax;
  };
  std::vector<element> elements;

  G4Cache<int> lastElement;

  // Nodes
  int nNodes;
  struct node {
   double x, y, z; // Coordinates
   double v; // Potential
   std::vector<double> w; // Weighting potentials
  };
  std::vector<node> nodes;

  // Materials
  int nMaterials;
  struct material {
    bool driftmedium;
    DriftMedium* medium;
  };
  std::vector<material> materials;

  // Bounding box
  bool hasBoundingBox;

  // Ranges
  double mapxmin, mapymin, mapzmin;
  double mapxmax, mapymax, mapzmax;
  double mapvmin, mapvmax;
  double regxmin, regymin, regzmin;
  double regxmax, regymax, regzmax;

  // Options
  bool checkMultipleElement; // Scan for multiple elements that contain a point
  bool warning; // Warnings flag
  bool debug;

  // Local coordinates
  // Calculate coordinates in linear tetrahedra
  int Coordinates12(double x, double y, double z, double& t1, double& t2,
                   double& t3, double& t4, int imap) const;
  // Calculate coordinates for curved quadratic tetrahedra
  int Coordinates13(double x, double y, double z, double& t1, double& t2,
                   double& t3, double& t4, double jac[4][4], double& det,
                   int imap) const;

  // Calculate Jacobian for curved quadratic tetrahedra
  void Jacobian13(int i, double t, double u, double v, double w, double& det,
                 double jac[4][4]) const;
  // Find the element for a point in curved quadratic tetrahedra
  int FindElement13(const double x, const double y, const double z, double& t1,
                   double& t2, double& t3, double& t4, double jac[4][4],
                   double& det);
  // Find the region for a given point.
  int FindRegion(const double x, const double y, const double z);

  static int ReadInteger(char* token, int def, bool& error);
  static double ReadDouble(char* token, double def, bool& error);

  // Calculate the bounding boxes of all elements after initialization
  void CalculateElementBoundingBoxes(void);
  // Calculate regions of elements after initialization
  void CalculateElementRegions(void);
};

#endif // FIELD_ELMER_MAP_H_
