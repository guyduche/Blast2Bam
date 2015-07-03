set term svg size 2048, 1280
set size 1, 1
set xlabel "position"
set ylabel "depth"
set style data lines
set output ARG1
plot @ARG2 using 2:3 title ARG3
