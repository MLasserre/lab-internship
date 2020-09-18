clear
reset

set terminal postscript eps enhanced "Helvetica" 20

unset title
unset key

set border linewidth 2.75
set xlabel font ",28"
set ylabel font ",26"

set autoscale
set xlabel 't'
set ylabel '<p^2>' 

plot 'qqvar5.dat' using 1:2 with lines linewidth 4

set output "qqvar5.eps"

replot
set terminal wxt enhanced
replot
