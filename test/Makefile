
.PHONY: clean destroy makeSrc
VPATH = ../bin/
FASTQ = ${HOME}/data/Pichon/180714_rand.fq.gz
DB = ${HOME}/data/Pichon/db/db.fasta


results.sam: makeSrc dbDict ${FASTQ} ${DB}
	gunzip -c ${FASTQ} | head -n 40000 | paste - - - - | cut -f 1,2 | cut -c2- |\
		sed 's/^/>/' | tr "\t" "\n" | blastn -db ${DB} -query - -outfmt 5 | ./../bin/prog - db.dict ${FASTQ} > $@
	rm -f db.dict

makeSrc:
	(cd ../bin; ${MAKE})

dbDict: ${DB}
	picard-tools CreateSequenceDictionary R=$< O=db.dict

resultsBWA.sam: ${FASTQ} ${DB}
	gunzip -c ${FASTQ} | head -n 40000 | bwa mem -p ${DB} - > $@
	
clean:
	rm -f *.sam db.dict
	
destroy: clean
	(cd ../bin; ${MAKE} clean)
