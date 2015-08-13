set term svg size 600, 400
set size 1, 1
set xlabel "word size"
set ylabel "number of reads mapped (%)"
unset key
set style data linespoint
set output ARG1
plot @ARG2 using @ARG3:((($3)*100)/@ARG4)
