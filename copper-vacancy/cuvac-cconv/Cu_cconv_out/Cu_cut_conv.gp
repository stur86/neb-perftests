
set xlabel 'Cut Off Energy (eV)'
set ylabel 'Deviation'
set title 'Convergence vs. plane wave cut off'
set xtics nomirror
set ytics nomirror

plot "Cu_cut_conv.dat" using 1:(($2-(-53786.80087575))*1.0) with linespoints pt 7 lc 1 ti "Final energy (eV)", "Cu_cut_conv.dat" using 1:(($3-(1.117226758227711))*1000.0) with linespoints pt 7 lc 2 ti "Max force (10^{-3} ev/Ang)", "Cu_cut_conv.dat" using 1:(($4-(0.0557769932244604))*100.0) with linespoints pt 7 lc 3 ti "Stress (norm, 10^{-2} GPa)", 0 lt 0 lc 0 notitle
pause -1 'Hit return to continue'
