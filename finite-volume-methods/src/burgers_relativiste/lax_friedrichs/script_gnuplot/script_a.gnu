set terminal gif animate

unset title
unset key

set border linewidth 2.75
set xlabel font ",22"
set ylabel font ",22"

set xlabel 'p'
set ylabel 'P(p,t)' 

set output "paquet.gif"
set xrange [-550:550]
set yrange [0:0.05]
i=0
imax=150
load 'script_b.gnu'
pause -1
replot
