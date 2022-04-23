/*  This is modified copy of Garfield++ TetrahedralTree with some optimizations.
 *  This helper class stores the mesh nodes and elements in an Octree data structure
 *  to optimize the element search operations
 */

#ifndef TetrahedralTree_H_
#define TetrahedralTree_H_

#include <vector>

#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <G4Cache.hh>

/* Brief Helper class for searches in field maps.
 * This class stores the mesh nodes and elements in an Octree data
 * structure to optimize the element search operations
 *
 * Author: Ali Sheharyar
 * Organization: Texas A&M University at Qatar
 */

class TetrahedralTree {
 public:
  TetrahedralTree(const G4ThreeVector& origin, const G4ThreeVector& halfDimension);
  ~TetrahedralTree();

  /// Insert a mesh node (a vertex/point) to the tree.
  void InsertMeshNode(const G4ThreeVector& point, const int index);
  /// Insert a mesh element with given bounding box and index to the tree.
  void InsertMeshElement(const double bb[6], const int index);
  /// Get all elements linked to a block corresponding to the given point.
  std::vector<int> GetElementsInBlock(const G4ThreeVector& point) const;
  /// Frees nodes which are unnecessary once all mesh elements are added to the tree.
  void Finalize(void);

 private:
  // Physical centre of this tree node.
  G4ThreeVector m_origin;
  // Half the width/height/depth of this tree node.
  G4ThreeVector m_halfDimension;
  // Storing min and max points for convenience
  G4ThreeVector m_min, m_max;

  // The tree has up to eight children and can additionally store
  // a list of mesh nodes and mesh elements.
  // Pointers to child octants.
  TetrahedralTree* children[8];

  // Children follow a predictable pattern to make accesses simple.
  // Here, - means less than 'origin' in that dimension, + means greater than.
  // child: 0 1 2 3 4 5 6 7
  // x:     - - - - + + + +
  // y:     - - + + - - + +
  // z:     - + - + - + - +

  std::vector<std::pair<G4ThreeVector, int> > nodes;
  std::vector<int> elements;

  static const size_t BlockCapacity = 10;

  // Check if the given box overlaps with this tree node.
  bool DoesBoxOverlap(const double bb[6]) const;

  int GetOctantContainingPoint(const G4ThreeVector& point) const;

  // Check if this tree node is a leaf or intermediate node.
  bool IsLeafNode() const;

  // Get a block containing the input point
  const TetrahedralTree* GetBlockFromPoint(const G4ThreeVector& point) const;

  // A helper function used by the function above.
  // Called recursively on the child nodes.
  const TetrahedralTree* GetBlockFromPointDescend(const G4ThreeVector& point) const;
};

#endif // TetrahedralTree_H_
