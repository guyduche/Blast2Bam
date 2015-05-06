
.PHONY: clean destroy makeSrc diff
VPATH = ../bin/
FASTQ = ${HOME}/data/HIV/test_1.gz ${HOME}/data/HIV/test_2.gz
DB = ${HOME}/data/HIV/db/hiv.fasta
SRA = ${HOME}/data/HIV/ERR656485.sra

diff: results.sam resultsBWA.sam
	$(foreach S,$^,cut -f 1,3,4 $S | grep -vF '@SQ' | grep -F 'gi|' | sed 's/	/~/' | LC_ALL=C sort -t '	' -k1,1 -k2,2  > $(addsuffix .tmp,$S); )
	join -t '	' -1 1 -2 1 $(addsuffix .tmp,$^) > same.read-chrom.txt
	join -t '	' -v1 -1 1 -2 1 $(addsuffix .tmp,$^) > only.prog.txt
	join -t '	' -v2 -1 1 -2 1 $(addsuffix .tmp,$^) > only.bwa.txt

# Prog
results.sam: ./../bin/prog ./../bin/fastqParser db.dict ${FASTQ} $(addsuffix .nin,${DB})
	./../bin/fastqParser ${FASTQ} | blastn -db ${DB} -query - -outfmt 5 | ./../bin/prog - db.dict ${FASTQ} > $@

resultsTest.sam: ./../bin/prog ./../bin/fastqParser db.dict ${FASTQ} $(addsuffix .nin,${DB})
	./../bin/fastqParser ${FASTQ} | blastn -db ${DB} -query - -outfmt 5 > test.xml

./../bin/prog ./../bin/fastqParser:
	(cd ../src; ${MAKE} )

db.dict: ${DB}
	picard-tools CreateSequenceDictionary R=$< O=$@

$(addsuffix .nin,${DB}) : ${DB}
	makeblastdb -in $< -dbtype nucl

# BWA
resultsBWA.sam: ${FASTQ} ${DB}
	bwa mem -t 5 ${DB} ${FASTQ} > $@

indexBWA: ${DB}
	bwa index ${DB}

# Fastq
fastq: ${SRA}
	(cd ${HOME}/data/HIV; fastq-dump --split-files --gzip ${SRA})

# Clean
clean:
	rm -f *.sam *.tmp *.txt db.dict
	
destroy: clean
	(cd ../src; ${MAKE} clean)
