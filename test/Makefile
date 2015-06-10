.PHONY: clean destroy bam_idx bamBWA_idx
.SHELL=/bin/bash
FASTQ = ${HOME}/data/nanopore/phage/ERR732557.fastq.gz
DB = ${HOME}/data/nanopore/phage/db/m13mp18.fasta
BWA = ${HOME}/programs/bwa-0.7.12/bwa
Blast = ${HOME}/programs/ncbi-blast-2.2.30+/bin/blastn
LAST = ${HOME}/programs/last-cea0cca54384

all: bam_idx bamBWA_idx

# Prog
bam_idx: resultsSorted.bam
	samtools index $<

resultsSorted.bam: results.sam
	samtools view -Sub $< | samtools sort - $(basename $@)

results.xml.gz: resultsSorted.sam

results.sam: ../bin/prog ../bin/fastqParser db.dict ${FASTQ} $(addsuffix .nin,${DB})
	../bin/fastqParser ${FASTQ} | ${Blast} -db ${DB} -num_threads 4 -word_size 8 -reward 1 -penalty -1 -gapopen 1 -gapextend 2 -outfmt 5 | \
	tee results.xml | ../bin/prog -W 40 -z -R '@RG	ID:foo	SM:sample' - db.dict ${FASTQ} > $@
	gzip --best -f results.xml

../bin/prog ../bin/fastqParser:
	(cd ../src; ${MAKE} )

db.dict: ${DB}
	picard-tools CreateSequenceDictionary R=$< O=$@

$(addsuffix .nin,${DB}): ${DB}
	makeblastdb -in $< -dbtype 'nucl'

# BWA
bamBWA_idx: resultsBWASorted.bam
	samtools index $<

resultsBWASorted.bam: resultsBWA.sam
	samtools view -Sub $< | samtools sort - $(basename $@)

resultsBWA.sam: ${FASTQ} $(addsuffix .bwt,${DB})
	${BWA} mem -t 4 -x 'ont2d' ${DB} ${FASTQ} > $@

$(addsuffix .bwt,${DB}): ${DB}
	${BWA} index ${DB}

# LAST
bamLAST_idx: resultsLASTSorted.bam
	samtools index $<

resultsLASTSorted.bam: resultsLAST.sam
	samtools view -Sub $< | samtools sort - $(basename $@)

resultsLAST.sam: lastdb ${FASTQ} db.dict
	gunzip -c ${FASTQ} | ${LAST}/src/lastal -Q 1 $< | ${LAST}/scripts/maf-convert -f db.dict sam > $@

lastdb: ${DB}
	${LAST}/src/lastdb $@ ${DB}

# Clean
clean:
	rm -f *.sam *.bam *.bai *.xml.gz db.dict lastdb*
	
destroy: clean
	(cd ../src; ${MAKE} clean)
