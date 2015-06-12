#! /bin/bash

for S in "$@" ; do
	samtools view $S | cut -f 1,3,4 | grep -vF '@SQ' | grep -F 'gi|' | sed 's/	/~/' | LC_ALL=C sort -t '	' -k 1,1 -k 2,2  > $S.tmp
done

join -t '	' -1 1 -2 1 $1.tmp $2.tmp > $(dirname "$1")/same.read-chrom.txt
join -t '	' -v 1 -1 1 -2 1 $1.tmp $2.tmp > $(dirname "$1")/only.prog.txt
join -t '	' -v 2 -1 1 -2 1 $1.tmp $2.tmp > $(dirname "$1")/only.bwa.txt
