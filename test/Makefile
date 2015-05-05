
.PHONY: clean destroy makeSrc
VPATH = ../bin/
FASTQ = ${HOME}/data/Ebola/SRR_1.gz ${HOME}/data/Ebola/SRR_2.gz
DB = ${HOME}/data/Ebola/db/ebolavirus.fasta

# Prog
results.sam: makeSrc dbDict ${FASTQ} ${DB}
	./../bin/fastqParser ${FASTQ} | blastn -db ${DB} -query - -outfmt 5 | ./../bin/prog - db.dict ${FASTQ} > $@
	rm -f db.dict

makeSrc:
	(cd ../src; ${MAKE})

dbDict: ${DB}
	picard-tools CreateSequenceDictionary R=$< O=db.dict

indexBlast: ${DB}
	(cd ${HOME}/data/Pichon/db; makeblastdb -in ${DB} -dbtype nucl)

# BWA
resultsBWA.sam: ${FASTQ} ${DB}
	bwa mem ${DB} ${FASTQ} > $@

indexBWA: ${DB}
	(cd ${HOME}/data/Pichon/db; bwa index ${DB})

# Clean
clean:
	rm -f *.sam db.dict
	
destroy: clean
	(cd ../src; ${MAKE} clean)
