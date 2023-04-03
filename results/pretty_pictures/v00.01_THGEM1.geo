// Gmsh project created on Thu Dec 01 2022
// All sizes are in mm
// THGEM #8 (CERN, 28% optical transparency).
pitch = 0.9;
rim = 0.1;
r_hole = 0.25; // Dielectric hole radius at the hole ends
r_center = 0.25; // Dielectric hole radius at the hole center
electrode = 0.03; // Electrode (copper) thiсkness
dielectric = 0.4; // Dielectric full thickness
cathode_z = 5; //from cathode to the bottom GEM electrode
anode_z = 5; //from cathode to the top GEM electrode

R = rim + r_hole; // Copper hole radius
TOP = cathode_z + electrode + dielectric/2;
BOT = -anode_z - electrode - dielectric/2;
x_size = pitch/4;
y_size = pitch/4 * Sqrt(3);
has_rim = Fabs(R - r_hole) > 1e-20;

// Anode points
Point(1) = {-x_size, -y_size, BOT, 1.0};
Point(2) = {-x_size, y_size, BOT, 1.0};
Point(3) = {x_size, y_size, BOT, 1.0};
Point(4) = {x_size, -y_size, BOT, 1.0};

// Bottom electrode's points
Point(5) = {-x_size, -y_size, -electrode - dielectric/2, 1.0};
Point(6) = {-x_size, y_size, -electrode - dielectric/2, 1.0};
Point(7) = {x_size, y_size, -electrode - dielectric/2, 1.0};
Point(8) = {x_size, -y_size, -electrode - dielectric/2, 1.0};

Point(9) = {-x_size + R, -y_size, -electrode - dielectric/2, 1.0};
Point(10) = {-x_size, -y_size + R, -electrode - dielectric/2, 1.0};
Point(11) = {x_size - R, y_size, -electrode - dielectric/2, 1.0};
Point(12) = {x_size, y_size - R, -electrode - dielectric/2, 1.0};

Point(13) = {-x_size, -y_size, -dielectric/2, 1.0};
Point(14) = {-x_size, y_size, -dielectric/2, 1.0};
Point(15) = {x_size, y_size, -dielectric/2, 1.0};
Point(16) = {x_size, -y_size, -dielectric/2, 1.0};

Point(17) = {-x_size + R, -y_size, -dielectric/2, 1.0};
Point(18) = {-x_size, -y_size + R, -dielectric/2, 1.0};
Point(19) = {x_size - R, y_size, -dielectric/2, 1.0};
Point(20) = {x_size, y_size - R, -dielectric/2, 1.0};

// Dielectric's points
Point(21) = {-x_size, -y_size, 0, 1.0};
Point(23) = {x_size, y_size, 0, 1.0};

Point(25) = {-x_size + r_center, -y_size, 0, 1.0};
Point(26) = {-x_size, -y_size + r_center, 0, 1.0};
Point(27) = {x_size - r_center, y_size, 0, 1.0};
Point(28) = {x_size, y_size - r_center, 0, 1.0};

If (has_rim)
  Point(49) = {-x_size + r_hole, -y_size, -dielectric/2, 1.0};
  Point(50) = {-x_size, -y_size + r_hole, -dielectric/2, 1.0};
  Point(51) = {x_size - r_hole, y_size, -dielectric/2, 1.0};
  Point(52) = {x_size, y_size - r_hole, -dielectric/2, 1.0};
  Point(53) = {-x_size + r_hole, -y_size, dielectric/2, 1.0};
  Point(54) = {-x_size, -y_size + r_hole, dielectric/2, 1.0};
  Point(55) = {x_size - r_hole, y_size, dielectric/2, 1.0};
  Point(56) = {x_size, y_size - r_hole, dielectric/2, 1.0};
EndIf

// Top electrode's points
Point(29) = {-x_size, -y_size, dielectric/2, 1.0};
Point(30) = {-x_size, y_size, dielectric/2, 1.0};
Point(31) = {x_size, y_size, dielectric/2, 1.0};
Point(32) = {x_size, -y_size, dielectric/2, 1.0};

Point(33) = {-x_size + R, -y_size, dielectric/2, 1.0};
Point(34) = {-x_size, -y_size + R, dielectric/2, 1.0};
Point(35) = {x_size - R, y_size, dielectric/2, 1.0};
Point(36) = {x_size, y_size - R, dielectric/2, 1.0};

Point(37) = {-x_size, -y_size, electrode + dielectric/2, 1.0};
Point(38) = {-x_size, y_size, electrode + dielectric/2, 1.0};
Point(39) = {x_size, y_size, electrode + dielectric/2, 1.0};
Point(40) = {x_size, -y_size, electrode + dielectric/2, 1.0};

Point(41) = {-x_size + R, -y_size, electrode + dielectric/2, 1.0};
Point(42) = {-x_size, -y_size + R, electrode + dielectric/2, 1.0};
Point(43) = {x_size - R, y_size, electrode + dielectric/2, 1.0};
Point(44) = {x_size, y_size - R, electrode + dielectric/2, 1.0};

// Cathode's points
Point(45) = {-x_size, -y_size, TOP, 1.0};
Point(46) = {-x_size, y_size, TOP, 1.0};
Point(47) = {x_size, y_size, TOP, 1.0};
Point(48) = {x_size, -y_size, TOP, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {6, 10};
Line(6) = {9, 8};
Line(7) = {8, 12};
Line(8) = {11, 6};
Circle(9) = {9, 5, 10};
Circle(10) = {12, 7, 11};
Line(11) = {16, 17};
Line(12) = {18, 14};
Line(13) = {14, 19};
Line(14) = {20, 16};
Circle(15) = {20, 15, 19};
Circle(16) = {17, 13, 18};
Line(17) = {8, 16};
Line(18) = {12, 20};
Line(19) = {11, 19};
Line(20) = {6, 14};
Line(21) = {10, 18};
Line(22) = {9, 17};
//+
Line(23) = {16, 32};
//Line(24) = {24, 32};
Line(25) = {14, 30};
//Line(26) = {22, 30};
Line(27) = {32, 33};
Line(28) = {30, 34};
Line(29) = {30, 35};
If (has_rim)
  Line(30) = {27, 51};
  Line(31) = {55, 27};
  Line(32) = {52, 28};
  Line(33) = {28, 56};
  Line(35) = {50, 26};
  Line(36) = {26, 54};
  Line(37) = {49, 25};
  Line(38) = {25, 53};
Else
  Line(30) = {27, 19};
  Line(31) = {35, 27};
  Line(32) = {20, 28};
  Line(33) = {28, 36};
  Line(35) = {18, 26};
  Line(36) = {26, 34};
  Line(37) = {17, 25};
  Line(38) = {25, 33};
EndIf
Line(34) = {36, 32};

Circle(39) = {26, 21, 25};
Circle(40) = {34, 29, 33};
Circle(41) = {28, 23, 27};
Circle(42) = {36, 31, 35};
//+
Line(43) = {44, 40};
Line(44) = {40, 41};
Line(45) = {42, 38};
Line(46) = {38, 43};
Line(47) = {35, 43};
Line(48) = {30, 38};
Line(49) = {34, 42};
Line(50) = {33, 41};
Line(51) = {32, 40};
Line(52) = {36, 44};
Circle(53) = {44, 39, 43};
Circle(54) = {41, 37, 42};
//+
Line(55) = {45, 46};
Line(56) = {46, 47};
Line(57) = {47, 48};
Line(58) = {48, 45};
//+
Line(59) = {45, 37};
Line(60) = {37, 29};
Line(61) = {29, 21};
Line(62) = {21, 13};
Line(63) = {13, 5};
Line(64) = {5, 1};
//+
Line(65) = {2, 6};
Line(66) = {38, 46};
//+
Line(67) = {47, 39};
Line(68) = {39, 31};
Line(69) = {31, 23};
Line(70) = {23, 15};
Line(71) = {15, 7};
Line(72) = {7, 3};
//+
Line(73) = {4, 8};
Line(74) = {40, 48};

If (has_rim)
  Circle(75) = {49, 13, 50};
  Circle(76) = {51, 15, 52};
  Circle(77) = {53, 29, 54};
  Circle(78) = {55, 31, 56};
  Line(79) = {49, 17};
  Line(80) = {50, 18};
  Line(81) = {51, 19};
  Line(82) = {52, 20};
  Line(83) = {53, 33};
  Line(84) = {54, 34};
  Line(85) = {55, 35};
  Line(86) = {56, 36};
EndIf

// Anode
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Bottom electrode
Line Loop(2) = {7, 18, 14, -17};
Plane Surface(2) = {2};
Line Loop(3) = {11, -22, 6, 17};
Plane Surface(3) = {3};
Line Loop(4) = {21, 12, -20, 5};
Plane Surface(4) = {4};
Line Loop(5) = {13, -19, 8, 20};
Plane Surface(5) = {5};
Line Loop(6) = {11, 16, 12, 13, -15, 14};
Plane Surface(6) = {6};
Line Loop(7) = {10, 8, 5, -9, 6, 7};
Plane Surface(7) = {7};
Line Loop(8) = {16, -21, -9, 22};
Surface(8) = {8};
Line Loop(9) = {15, -19, -10, 18};
Surface(9) = {9};

//dielectic
If (has_rim)
  // Dielectric vertical sides
  Line Loop(10) = {33, 86, 34, -23, -14, -82, 32};
  Line Loop(11) = {11, -79, 37, 38, 83, -27, -23};
  Line Loop(12) = {28, -84, -36, -35, 80, 12, 25};
  Line Loop(13) = {29, -85, 31, 30, 81, -13, 25};
  // Dielectric hole 
  Line Loop(15) = {32, 41, 30, 76};
  Line Loop(16) = {33, -78, 31, -41};
  Line Loop(17) = {37, -39, -35, -75};
  Line Loop(18) = {38, 77, -36, 39};
Else
  // Dielectric vertical sides
  Line Loop(10) = {33, 34, -23, -14, 32};
  Line Loop(11) = {11, 37, 38, -27, -23};
  Line Loop(12) = {28, -36, -35, 12, 25};
  Line Loop(13) = {29, 31, 30, -13, 25};
  // Dielectric hole
  Line Loop(15) = {32, 41, 30, -15};
  Line Loop(16) = {33, 42, 31, -41};
  Line Loop(17) = {37, -39, -35, -16};
  Line Loop(18) = {38, -40, -36, 39};
EndIf
Plane Surface(10) = {10};
Plane Surface(11) = {11};
Plane Surface(12) = {12};
Plane Surface(13) = {13};

Surface(15) = {15};
Surface(16) = {16};
Surface(17) = {17};
Surface(18) = {18};

Line Loop(14) = {27, -40, -28, 29, -42, 34};
Plane Surface(14) = {14};

If (has_rim)
  Line Loop(31) = {-42, -86, -78, 85};
  Plane Surface(31) = {31};
  Line Loop(32) = {-15, -82, -76, 81};
  Plane Surface(32) = {32};
  Line Loop(33) = {-40, -84, -77, 83};
  Plane Surface(33) = {33};
  Line Loop(34) = {16, -80, -75, 79};
  Plane Surface(34) = {34};
EndIf

// Upper electrode
Line Loop(19) = {44, -50, -27, 51};
Plane Surface(19) = {19};
Line Loop(20) = {45, -48, 28, 49};
Plane Surface(20) = {20};
Line Loop(21) = {46, -47, -29, 48};
Plane Surface(21) = {21};
Line Loop(22) = {43, -51, -34, 52};
Plane Surface(22) = {22};
Line Loop(23) = {44, 54, 45, 46, -53, 43};
Plane Surface(23) = {23};
Line Loop(24) = {53, -47, -42, 52};
Surface(24) = {24};
Line Loop(25) = {54, -49, 40, 50};
Surface(25) = {25};

// Cathode
Line Loop(26) = {58, 55, 56, 57};
Plane Surface(26) = {26};

// Gas
If (has_rim)
  Line Loop(27) = {2, -72, -71, -70, -69, -68, -67, -56, -66, 46, -47, -85, 31, 30, 81, -19, 8, -65};
  Line Loop(28) = {58, 59, 60, 61, 62, 63, 64, -4, 73, -6, 22, -79, 37, 38, 83, 50, -44, 74};
  Line Loop(29) = {64, 1, 65, 5, 21, -80, 35, 36, 84, 49, 45, 66, -55, 59, 60, 61, 62, 63};
  Line Loop(30) = {67, 68, 69, 70, 71, 72, 3, 73, 7, 18, -82, 32, 33, 86, 52, 43, 74, -57};
Else
  Line Loop(27) = {2, -72, -71, -70, -69, -68, -67, -56, -66, 46, -47, 31, 30, -19, 8, -65};
  Line Loop(28) = {58, 59, 60, 61, 62, 63, 64, -4, 73, -6, 22, 37, 38, 50, -44, 74};
  Line Loop(29) = {64, 1, 65, 5, 21, 35, 36, 49, 45, 66, -55, 59, 60, 61, 62, 63};
  Line Loop(30) = {67, 68, 69, 70, 71, 72, 3, 73, 7, 18, 32, 33, 52, 43, 74, -57};
EndIf
Plane Surface(27) = {27};
Plane Surface(28) = {28};
Plane Surface(29) = {29};
Plane Surface(30) = {30};

// Gas
If (has_rim)
  Surface Loop(1) = {27, 1, 29, 28, 26, 30, 7, 9, 8, 17, 18, 23, 25, 24, 16, 15, 31, 32, 33, 34};
Else
  Surface Loop(1) = {27, 1, 29, 28, 26, 30, 7, 9, 8, 17, 18, 23, 25, 24, 16, 15};
EndIf
Volume(1) = {1};

// Dielectric
If (has_rim)
  Surface Loop(2) = {15, 16, 13, 12, 6, 14, 18, 17, 11, 10, 31, 32, 33, 34};
Else
  Surface Loop(2) = {15, 16, 13, 12, 6, 14, 18, 17, 11, 10};
EndIf
Volume(2) = {2};
// Bottom electrode
Surface Loop(3) = {7, 9, 8, 6, 5, 4, 3, 2};
Volume(3) = {3};
// Top electrode
Surface Loop(4) = {25, 23, 24, 20, 21, 22, 19, 14};
Volume(4) = {4};

Field[1] = Box;
Field[1].VIn = 0.08;
Field[1].VOut = 10;
Field[1].XMax = 1;
Field[1].XMin = -1;
Field[1].YMax = 1;
Field[1].YMin = -1;
Field[1].ZMax = 100;
Field[1].ZMin = -100;
//+
Field[2] = Box;
Field[2].VIn = 0.015;
Field[2].VOut = 10;
Field[2].XMax = 1;
Field[2].XMin = -1;
Field[2].YMax = 1;
Field[2].YMin = -1;
Field[2].ZMax = 0.3;
Field[2].ZMin = -0.3;
//+
Field[3] = Box;
Field[3].VIn = 0.03;
Field[3].VOut = 10;
Field[3].XMax = 1;
Field[3].XMin = -1;
Field[3].YMax = 1;
Field[3].YMin = -1;
Field[3].ZMax = 0.5;
Field[3].ZMin = -0.5;

Field[4] = Min;
Field[4].FieldsList = {1, 2, 3};
Background Field = 4;

Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.ElementOrder = 2;
Mesh.HighOrderFixBoundaryNodes = 1;
Mesh.HighOrderOptimize = 2;
Mesh.HighOrderThresholdMin = 0.1; // default
Mesh.HighOrderThresholdMax = 2.0; // default
Mesh.OptimizeThreshold = 0.1;
Mesh.SaveAll = 1;
