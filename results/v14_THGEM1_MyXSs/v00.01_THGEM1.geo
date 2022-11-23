// Gmsh project created on Wed Nov 10 11:47:30 2021
//All sizes are in mm
cathode_z = 5; //from cathode to bottom THGEM electrode
anode_z = 5; //from cathode to top THGEM electrode
dielec = 0.4; //total dielectric thickness
electrode = 0.03; //electrode thickness
rim = 0.1;
hole_r = 0.25;
pitch = 0.9;

TOP = cathode_z + electrode + dielec/2;
BOT = -anode_z - electrode - dielec/2;
X_SZ = pitch/4;
Y_SZ = pitch/4 * Sqrt(3);
R = rim + hole_r;
r = hole_r;

Point(1) = {X_SZ, Y_SZ, TOP, 1.0};
Point(2) = {X_SZ, Y_SZ, BOT, 1.0};
Point(3) = {-X_SZ, Y_SZ, TOP, 1.0};
Point(4) = {-X_SZ, Y_SZ, BOT, 1.0};
Point(5) = {-X_SZ, -Y_SZ, BOT, 1.0};
Point(6) = {-X_SZ, -Y_SZ, TOP, 1.0};
Point(7) = {X_SZ, -Y_SZ, TOP, 1.0};
Point(8) = {X_SZ, -Y_SZ, BOT, 1.0};
//+
Point(9) = {X_SZ, Y_SZ, electrode + dielec/2, 1.0};
Point(10) = {X_SZ, Y_SZ, dielec/2, 1.0};
Point(11) = {X_SZ, Y_SZ, -dielec/2, 1.0};
Point(12) = {X_SZ, Y_SZ, -electrode - dielec/2, 1.0};
Point(13) = {-X_SZ, Y_SZ, -electrode - dielec/2, 1.0};
Point(14) = {-X_SZ, Y_SZ, -dielec/2, 1.0};
Point(15) = {-X_SZ, Y_SZ, dielec/2, 1.0};
Point(16) = {-X_SZ, Y_SZ, electrode + dielec/2, 1.0};
//+
Point(17) = {-X_SZ, -Y_SZ, electrode + dielec/2, 1.0};
Point(18) = {-X_SZ, -Y_SZ, dielec/2, 1.0};
Point(19) = {-X_SZ, -Y_SZ, -dielec/2, 1.0};
Point(20) = {-X_SZ, -Y_SZ, -electrode - dielec/2, 1.0};
Point(21) = {X_SZ, -Y_SZ, -electrode - dielec/2, 1.0};
Point(22) = {X_SZ, -Y_SZ, -dielec/2, 1.0};
Point(23) = {X_SZ, -Y_SZ, dielec/2, 1.0};
Point(24) = {X_SZ, -Y_SZ, electrode + dielec/2, 1.0};
//+
Point(25) = {X_SZ - R, Y_SZ, electrode + dielec/2, 1.0};
Point(26) = {X_SZ - R, Y_SZ, dielec/2, 1.0};
Point(27) = {X_SZ - r, Y_SZ, dielec/2, 1.0};
Point(28) = {X_SZ - r, Y_SZ, -dielec/2, 1.0};
Point(29) = {X_SZ - R, Y_SZ, -dielec/2, 1.0};
Point(30) = {X_SZ - R, Y_SZ, -electrode - dielec/2, 1.0};
Point(31) = {X_SZ, Y_SZ - R, electrode + dielec/2, 1.0};
Point(32) = {X_SZ, Y_SZ - R, dielec/2, 1.0};
Point(33) = {X_SZ, Y_SZ - r, dielec/2, 1.0};
Point(34) = {X_SZ, Y_SZ - r, -dielec/2, 1.0};
Point(35) = {X_SZ, Y_SZ - R, -dielec/2, 1.0};
Point(36) = {X_SZ, Y_SZ - R, -electrode - dielec/2, 1.0};
//+
Point(37) = {-X_SZ + R, -Y_SZ, electrode + dielec/2, 1.0};
Point(38) = {-X_SZ + R, -Y_SZ, dielec/2, 1.0};
Point(39) = {-X_SZ + r, -Y_SZ, dielec/2, 1.0};
Point(40) = {-X_SZ + r, -Y_SZ, -dielec/2, 1.0};
Point(41) = {-X_SZ + R, -Y_SZ, -dielec/2, 1.0};
Point(42) = {-X_SZ + R, -Y_SZ, -electrode - dielec/2, 1.0};
Point(43) = {-X_SZ, -Y_SZ + R, -electrode - dielec/2, 1.0};
Point(44) = {-X_SZ, -Y_SZ + R, -dielec/2, 1.0};
Point(45) = {-X_SZ, -Y_SZ + r, -dielec/2, 1.0};
Point(46) = {-X_SZ, -Y_SZ + r, dielec/2, 1.0};
Point(47) = {-X_SZ, -Y_SZ + R, dielec/2, 1.0};
Point(48) = {-X_SZ, -Y_SZ + R, electrode + dielec/2, 1.0};



//+
Line(1) = {7, 6};
//+
Line(2) = {6, 3};
//+
Line(3) = {3, 1};
//+
Line(4) = {1, 7};
//+
Line(5) = {1, 2};
//+
Line(6) = {5, 6};
//+
Line(7) = {3, 16};
//+
Line(8) = {16, 15};
//+
Line(9) = {15, 14};
//+
Line(10) = {14, 13};
//+
Line(11) = {13, 4};
//+
Line(12) = {8, 21};
//+
Line(13) = {21, 22};
//+
Line(14) = {22, 23};
//+
Line(15) = {23, 24};
//+
Line(16) = {24, 7};
//+
Line(17) = {8, 5};
//+
Line(18) = {5, 4};
//+
Line(19) = {4, 2};
//+
Line(20) = {2, 8};
//+
Line(21) = {16, 25};
//+
Line(22) = {25, 26};
//+
Line(23) = {26, 15};
//+
Line(24) = {16, 48};
//+
Line(25) = {48, 47};
//+
Line(26) = {15, 47};
//+
Line(27) = {47, 46};
//+
Line(28) = {29, 14};
//+
Line(29) = {29, 30};
//+
Line(30) = {30, 13};
//+
Line(31) = {14, 44};
//+
Line(32) = {44, 43};
//+
Line(33) = {43, 13};
//+
Line(34) = {44, 45};
//+
Line(35) = {28, 29};
//+
Line(36) = {26, 27};
//+
Line(37) = {27, 28};
//+
Line(38) = {46, 45};
//+
Line(39) = {24, 37};
//+
Line(40) = {37, 38};
//+
Line(41) = {38, 23};
//+
Line(42) = {24, 31};
//+
Line(43) = {31, 32};
//+
Line(44) = {32, 23};
//+
Line(45) = {32, 33};
//+
Line(46) = {38, 39};
//+
Line(47) = {41, 22};
//+
Line(48) = {41, 42};
//+
Line(49) = {42, 21};
//+
Line(50) = {21, 36};
//+
Line(51) = {36, 35};
//+
Line(52) = {35, 34};
//+
Line(53) = {35, 22};
//+
Line(54) = {41, 40};
//+
Line(55) = {40, 39};
//+
Line(56) = {33, 34};
//+
Circle(57) = {28, 11, 34};
//+
Circle(58) = {30, 12, 36};
//+
Circle(59) = {27, 10, 33};
//+
Circle(60) = {25, 9, 31};
//+
Circle(61) = {46, 18, 39};
//+
Circle(62) = {48, 17, 37};
//+
Circle(63) = {45, 19, 40};
//+
Circle(64) = {43, 20, 42};
//+
Circle(65) = {26, 10, 32};
//+
Circle(66) = {47, 18, 38};
//+
Circle(67) = {41, 19, 44};
//+
Circle(68) = {29, 11, 35};
//+
Line Loop(1) = {18, 19, 20, 17};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {1, 2, 3, 4};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {4, -16, 42, 43, 45, 56, -52, -51, -50, -12, -20, -5};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {18, -11, -33, -32, 34, -38, -27, -25, -24, -7, -2, -6};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {1, -6, -17, 12, -49, -48, 54, 55, -46, -40, -39, 16};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {19, -5, -3, 7, 21, 22, 36, 37, 35, 29, 30, 11};
//+
Plane Surface(6) = {6};
//+
Line Loop(7) = {45, 56, -52, 53, 14, -44};
//+
Plane Surface(7) = {7};
//+
Line Loop(8) = {27, 38, -34, -31, -9, 26};
//+
Plane Surface(8) = {8};
//+
Line Loop(9) = {55, -46, 41, -14, -47, 54};
//+
Plane Surface(9) = {9};
//+
Line Loop(10) = {36, 37, 35, 28, -9, -23};
//+
Plane Surface(10) = {10};
//+
Line Loop(11) = {60, -42, 39, -62, -24, 21};
//+
Plane Surface(11) = {11};
//+
Line Loop(12) = {21, 22, 23, -8};
//+
Plane Surface(12) = {12};
//+
Line Loop(13) = {26, -25, -24, 8};
//+
Plane Surface(13) = {13};
//+
Line Loop(14) = {40, 41, 15, 39};
//+
Plane Surface(14) = {14};
//+
Line Loop(15) = {44, 15, 42, 43};
//+
Plane Surface(15) = {15};
//+
Line Loop(16) = {66, 41, -44, -65, 23, 26};
//+
Plane Surface(16) = {16};
//+
Line Loop(17) = {66, -40, -62, 25};
//+
Surface(17) = {17};
//+
Line Loop(18) = {60, 43, -65, -22};
//+
Surface(18) = {18};
//+
Line Loop(19) = {68, 53, -47, 67, -31, -28};
//+
Plane Surface(19) = {19};
//+
Line Loop(20) = {33, -30, 58, -50, -49, -64};
//+
Plane Surface(20) = {20};
//+
Line Loop(21) = {28, 10, -30, -29};
//+
Plane Surface(21) = {21};
//+
Line Loop(22) = {33, -10, 31, 32};
//+
Plane Surface(22) = {22};
//+
Line Loop(23) = {47, -13, -49, -48};
//+
Plane Surface(23) = {23};
//+
Line Loop(24) = {50, 51, 53, -13};
//+
Plane Surface(24) = {24};
//+
Line Loop(25) = {67, 32, 64, -48};
//+
Surface(25) = {25};
//+
Line Loop(26) = {68, -51, -58, -29};
//+
Surface(26) = {26};
//+
Line Loop(27) = {34, 63, -54, 67};
//+
Plane Surface(27) = {27};
//+
Line Loop(28) = {57, -52, -68, -35};
//+
Plane Surface(28) = {28};
//+
Line Loop(29) = {36, 59, -45, -65};
//+
Plane Surface(29) = {29};
//+
Line Loop(30) = {61, -46, -66, 27};
//+
Plane Surface(30) = {30};
//+
Line Loop(31) = {61, -55, -63, -38};
//+
Surface(31) = {31};
//+
Line Loop(32) = {59, 56, -57, -37};
//+
Surface(32) = {32};
//+
Surface Loop(1) = {5, 2, 4, 1, 6, 3, 31, 30, 27, 29, 32, 28, 17, 11, 18, 25, 20, 26};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {9, 7, 31, 30, 27, 32, 29, 28, 10, 8, 16, 19};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {11, 18, 17, 15, 14, 16, 13, 12};
//+
Volume(3) = {3};
//+
Surface Loop(4) = {22, 21, 26, 20, 25, 23, 24, 19};
//+
Volume(4) = {4};

Field[1] = Box;
Field[1].VIn = 0.10;
Field[1].VOut = 10;
Field[1].XMax = 1;
Field[1].XMin = -1;
Field[1].YMax = 1;
Field[1].YMin = -1;
Field[1].ZMax = 50;
Field[1].ZMin = -50;
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
