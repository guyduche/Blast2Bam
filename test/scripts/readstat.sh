#! /bin/bash

for S in "$@" ; do
	gunzip -c $S | paste - - - - | cut -d "	" -f 2 > "$S".tmp
	nbReads=$(wc -l "$S".tmp | cut -d " " -f 1)
	nbBases=$(tr -d "\n" < "$S".tmp | wc -m)
	echo "$S :"
	echo "Number of bases : $nbBases"
	echo "Number of reads : $nbReads"
	echo "Average read size : $((nbBases / nbReads))"
	rm -f "$S".tmp
done
