
set xlabel 'Cut Off Energy (eV)'
set ylabel 'Deviation'
set title 'Convergence vs. plane wave cut off'
set xtics nomirror
set ytics nomirror

plot "AlAu_cut_conv.dat" using 1:(($2-(-15456.10149622))*1.0) with linespoints pt 7 lc 1 ti "Final energy (eV)", "AlAu_cut_conv.dat" using 1:(($3-(0.9089614798218899))*10.0) with linespoints pt 7 lc 2 ti "Max force (10^{-1} ev/Ang)", "AlAu_cut_conv.dat" using 1:(($4-(0.014784595462950486))*100.0) with linespoints pt 7 lc 3 ti "Stress (norm, 10^{-2} GPa)", 0 lt 0 lc 0 notitle
pause -1 'Hit return to continue'
