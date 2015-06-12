#! /bin/bash

for S in "$@" ; do
	gunzip -c $S | grep -vF "@ERR" | paste - - - | cut -d "	" -f 1 > "$S".tmp
	nbReads=$(wc -l "$S".tmp | cut -d " " -f 1)
	nbBases=$(tr -d "\n" < "$S".tmp | wc -m)
	echo "$S :"
	echo "Number of bases : $nbBases"
	echo "Number of reads : $nbReads"
	echo "Average read size : $((nbBases / nbReads))"
	rm "$S".tmp
done
