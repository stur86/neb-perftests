
set xlabel 'K-Points Grid'
set ylabel 'Deviation'
set title 'Convergence vs. k-points'
set xtics nomirror
set ytics nomirror
set xtics ("2x2x1" 4, "4x4x2" 32, "7x7x3" 147) rotate by 45 right
plot "AlAu_kpn_conv.dat" using 1:(($2-(-15448.11559333))*1.0) with linespoints pt 7 lc 1 ti "Final energy (eV)", "AlAu_kpn_conv.dat" using 1:(($3-(0.7710368080059473))*1.0) with linespoints pt 7 lc 2 ti "Max force (10^{0} ev/Ang)", "AlAu_kpn_conv.dat" using 1:(($4-(0.0230362759745022))*10.0) with linespoints pt 7 lc 3 ti "Stress (norm, 10^{-1} GPa)", 0 lt 0 lc 0 notitle
pause -1 'Hit return to continue'
