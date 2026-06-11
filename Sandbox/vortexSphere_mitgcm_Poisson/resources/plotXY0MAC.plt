set xlabel "SEDxLabel"
set ylabel "SEDyLabel"
m="../resources/SEDfileName.txt"
#set terminal x11 0
set terminal png
set output '../results/SEDfileName.png'
#set nokey
set grid
#set title 'SEDPlotTitle'
plot m using 1:2 with linespoints title 'SEDyLabel'
