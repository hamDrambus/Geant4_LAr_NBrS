set xlabel "e path [mm]" 
set xrange [0:10]
set grid

prefix = "../tests/test_10_e_tracks/e_track_T0_E"
suffix = ".txt"

#set ylabel "e z position [mm]"
#set yrange [-3.5:0.75]
#plot for [k=0:19] sprintf("%s%i%s", prefix, k, suffix) u 1:4 title sprintf("%i%s", k, suffix)

#Field amplitude vs length
set ylabel "Electric field [kV/cm]"
set yrange [0:100]
plot for [k=0:19] sprintf("%s%i%s", prefix, k, suffix) u 1:8 title sprintf("%i%s", k, suffix)

#set ylabel "Drift time [us]"
#set yrange [0:0.5]
#plot for [k=0:19] sprintf("%s%i%s", prefix, k, suffix) u 1:9 title sprintf("%i%s", k, suffix)

#set ylabel "Electric field Z [kV/cm]"
#set yrange [-100:100]
#plot for [k=0:19] sprintf("%s%i%s", prefix, k, suffix) u 1:7 title sprintf("%i%s", k, suffix)

#set xlabel "e z position [mm]"
#set xrange [-3:1]
#set ylabel "Electric field Z [kV/cm]"
#set yrange [-100:100]
#plot for [k=0:19] sprintf("%s%i%s", prefix, k, suffix) u 4:7 title sprintf("%i%s", k, suffix)

#set xlabel "e z position [mm]"
#set xrange [-3:1]
#set ylabel "Electric field [kV/cm]"
#set yrange [0:100]
#plot for [k=0:19] sprintf("%s%i%s", prefix, k, suffix) u 4:8 title sprintf("%i%s", k, suffix)

#set xlabel "e z position [mm]"
#set xrange [-3:1]
#set ylabel "e x position [mm]"
#set yrange [-1:1]
#plot for [k=0:19] sprintf("%s%i%s", prefix, k, suffix) u 4:2 title sprintf("%i%s", k, suffix)
