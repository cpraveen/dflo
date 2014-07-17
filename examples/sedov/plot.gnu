set xran[0:1.1]
set ylabel 'Density'
set xlabel 'x'
set key left top
p 'spherical_standard_omega0p00.dat' u 2:3 t 'Exact' w l lw 2, \
  'Density.curve' u ($1-2):2 t 'dflo' w p lw 2
set term postscript enhanced
set out 'sedov.eps'
replot
