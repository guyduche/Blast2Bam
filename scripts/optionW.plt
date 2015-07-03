set term svg size 600, 400
set size 1, 1
set xlabel "minimum alignment length"
set ylabel "number of reads mapped"
unset key
set style data linespoint
set output ARG1
plot @ARG2 using ($0 * @ARG3):3
