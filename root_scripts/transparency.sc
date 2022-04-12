# This is gnuplot script
# Geant4 resulting data must be present in separate file

# THGEM dielectric thickness [mm]
d = 0.4
# THGEM hole radius [mm]
r = 0.25
th(x) = x * pi / 180
# distance between circles
D(x) = tan(th(x)) * d
# Area of circle intersection
S(x) = abs(D(x)) > 2*r ? 0 : (2*r*r*acos(D(x)/(2*r)) - 0.5 * D(x) * sqrt(4*r*r - D(x)*D(x)))

# THGEM full thickness [mm]
d1 = 0.46
# THGEM hole+rim radius [mm]
r1 = 0.35
# distance between circles
D1(x) = tan(th(x)) * d1
# Area of circle intersection
S1(x) = abs(D1(x)) > 2*r1 ? 0 : (2*r1*r1*acos(D1(x)/(2*r1)) - 0.5 * D1(x) * sqrt(4*r1*r1 - D1(x)*D1(x)))

# THGEM theoretical transparency as a funciton of angle in degrees
# Both limitations by dielectric and copper are considered
Fr(x) = S(x) / (pi*r*r)
Fr1(x) = S1(x) / (pi*r1*r1)
Tr(x) = 28.00 *(Fr(x) > Fr1(x) ? Fr1(x) : Fr(x))

set xlabel "Angle [degree]"
set ylabel "THGEM1 transparency [%]"
set xrange [0:90]
plot Tr(x) title "Theoretical THGEM transparency" with lines
replot "../tests/test_06_SiPM_THGEM1_shading/results.txt" u 1:(100*$2):(100*$3) with errorbars pt 5 title "MC results"
