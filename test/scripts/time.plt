set term svg size 300, 200
set size 1, 1
set xlabel "BLAST seed length"
set ylabel "Time (s)"
unset key
set style data linespoint
set output ARG1
plot @ARG2 using @ARG3:1
