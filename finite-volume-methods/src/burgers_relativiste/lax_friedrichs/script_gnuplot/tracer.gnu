set xrange[0:300]
set yrange[0:250]
plot 'Var5.dat' using 1:2 pointsize 0.1
set terminal png  
set output "resultat.png"
replot
set terminal x11
