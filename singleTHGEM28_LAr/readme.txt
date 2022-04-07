version legend:
v[geometry version].[meshing version].[elmer solution version]
v00.x.x -  0.5mm thickness 28% THGEM0, dl from THGEM to cathode/anonde == 3mm
vx.01.x - normal meshing, finer mesh near TGHEM0 hole
Electric fields and resistances are according to 2022.02.03 measurements with THGEM1 in LAr and
pulsed X-ray source
vx.x.01 - V0 = 20kV (drift field 6.81 kV/cm), V1(on THGEM1 divider) = 6180 V
vx.x.02 - V0 = 20kV (drift field 6.81 kV/cm), V1(on THGEM1 divider) = 5993 V
vx.x.03 - V0 = 20kV (drift field 6.81 kV/cm), V1(on THGEM1 divider) = 5738 V
vx.x.04 - V0 = 20kV (drift field 6.81 kV/cm), V1(on THGEM1 divider) = 5297 V
vx.x.05 - V0 = 20kV (drift field 6.81 kV/cm), V1(on THGEM1 divider) = 4856 V
vx.x.06 - V0 = 20kV (drift field 6.81 kV/cm), V1(on THGEM1 divider) = 4413 V
vx.x.07 - V0 = 20kV (drift field 6.81 kV/cm), V1(on THGEM1 divider) = 3972 V
vx.x.08 - V0 = 20kV (drift field 6.81 kV/cm), V1(on THGEM1 divider) = 3531 V
vx.x.09 - V0 = 20kV (drift field 6.81 kV/cm), V1(on THGEM1 divider) = 3090 V
vx.x.10 - V0 = 20kV (drift field 6.81 kV/cm), V1(on THGEM1 divider) = 2648 V
vx.x.11 - V0 = 20kV (drift field 6.81 kV/cm), V1(on THGEM1 divider) = 2206 V
vx.x.11 - V0 = 20kV (drift field 6.81 kV/cm), V1(on THGEM1 divider) = 1765 V

v00.00:
==================================================================================================
+-----------+        .  .      +------+.......+     ----.
|           |          /|\     |------|-------|  .  ____| t
+___.        \          |      |      |-    - | / \ 
     \        \         |a     |      |_  -  _|  T
      \         *____+  |  X   |______|___-___| \./_____.
       *             |  |      +______+.......+    _____| t
       |             |  |
.      +_____________+ \./
<--------a┐/3'--------->
			Y
28% Optical transparency
r = 0.25 mm
R = 0.35 mm
a = P/2 = 0.45 mm
t = 0.03 mm
T = 0.4 mm - FR4
dz_bot = 3 mm - the distance from the bottom THGEM1 surface to the bottom plane creating the external field (cathode)
dz_top = 3 mm - the distance from the top THGEM1 surface to the top plane creating the external field (anode)
==================================================================================================
Surfaces:
Cathode (Z-): 1
Anode (Z+): 2
X+ LAr: 3	X+ diel: 7
X- LAr: 4	X- diel: 8
Y- LAr: 5	Y- diel: 9
Y+ LAr: 6	Y+ diel: 10
Cu top: 11, 12, 13, 14, 15, 16, 17, 18
Cu bot: 19, 20, 21, 22, 23, 24, 25, 26

Volumes:	1 - LAr	
		2 - dielectric				
		3 - top Cu
		4 - bot Cu
e_gas_Ar = 1.0
e_LAr = 1.54
e_FR4 (fiber glass) = 4.4
e_acrylic = 3.6

---------------------			0 kV	 Ground
		          		      ______		.
						 ____|__R4__|__0 V1 .
--------------------- <-.|.		1.2 MOhm	.
--||-----||-----||---  	|_|R5 = 8.6 MOhm	. THGEM1
--------------------- <- |					.
						 |   ______			.
						 |__|__R6__|__0  Ground
			     			1 MOhm			.
											.
		          		      ______		.
_____________________    ____|__R4__|__0 Ground
~	~	~  	~   ~   ~	 |	600 MOhm		.
--------------------- <-.|.					.
--||-----||-----||---  	|_|R3 = 0 Ohm(mesh) . THGEM0 or mesh
--------------------- <- |					.
~     ~		~			.|.					.
	~	~				| |R2 = 120 MOhm	.
   ~	     ~     ~	|_|  ______			.
--------------------- <--|__|__R1__|__0  V0	.
			     			80 MOhm			.


==================================================================================================
==================================================================================================
==================================================================================================
v00.01 - same as v00.00, but with meshing parameters:
Field[1] = Box;
Field[1].VIn = 0.1;
Field[1].VOut = 10;
Field[1].XMax = 1;
Field[1].XMin = -1;
Field[1].YMax = 1;
Field[1].YMin = -1;
Field[1].ZMax = 5;
Field[1].ZMin = -5;
//+
Field[2] = Box;
Field[2].VIn = 0.01;
Field[2].VOut = 10;
Field[2].XMax = 1;
Field[2].XMin = -1;
Field[2].YMax = 1;
Field[2].YMin = -1;
Field[2].ZMax = 0.5;
Field[2].ZMin = -0.5;
//+
Field[3] = Min;
Field[3].FieldsList = {1, 2};

Background Field = 3;

Elmer v01:
	V0 = 20kV
	V1 = 6180 V  Vthgem1 = 6180*8.6/10.8 = 4921 V
	V_bot_THGEM1 = 6180 * 1/10.8 = +572.2 V
	V_top_THGEM1 = 6180 * 9.6/10.8 = +5493.3 V
	Edrift = (20 kV * 600/800 + 572.2 V) /2.2 =  7.08 kV/cm = 708 V/mm
	Vcathode = 708 * 3mm = -2124 V;
	Einduction = 5493.3 V / 0.5 cm = 10.99 kV/cm = 1099 V/mm
	V_anode = 5493.3 V - 1099 * 3 = +2196.3 V


