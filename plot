set title "Spline on large data"
plot [-5:17] [-10000:60000] "large_data.dat"  with points
replot "pts" with lines 
set output "plot.ps"
set terminal postscript enhanced color landscape
replot
