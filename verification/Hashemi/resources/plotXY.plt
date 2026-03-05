set xlabel "theta"
set ylabel "Coriolis Term"
m="../resources/coriolisTerm.txt"
#set terminal x11 0
set terminal png
set output '../results/coriolisTerm.png'
#set nokey
set grid
#set title 'SEDPlotTitle'
plot m using 1:2 with linespoints title 'Coriolis Term'
