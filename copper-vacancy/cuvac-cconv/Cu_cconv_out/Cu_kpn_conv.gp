
set xlabel 'K-Points Grid'
set ylabel 'Deviation'
set title 'Convergence vs. k-points'
set xtics nomirror
set ytics nomirror
set xtics ("2x2x2" 8, "3x3x3" 27, "4x4x4" 64, "5x5x5" 125) rotate by 45 right
plot "Cu_kpn_conv.dat" using 1:(($2-(-53787.70653728))*1.0) with linespoints pt 7 lc 1 ti "Final energy (eV)", "Cu_kpn_conv.dat" using 1:(($3-(1.1612582919402556))*10.0) with linespoints pt 7 lc 2 ti "Max force (10^{-1} ev/Ang)", "Cu_kpn_conv.dat" using 1:(($4-(0.04548436744022302))*100.0) with linespoints pt 7 lc 3 ti "Stress (norm, 10^{-2} GPa)", 0 lt 0 lc 0 notitle
pause -1 'Hit return to continue'
