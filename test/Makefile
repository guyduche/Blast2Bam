.PHONY: clean destroy bam_idx bamBWA_idx
.SHELL=/bin/bash
BINDIR = ../bin/
FASTQ = ${HOME}/data/nanopore/phage/ERR732557.fastq.gz
DB = ${HOME}/data/nanopore/phage/db/m13mp18.fasta
BWA = ${HOME}/programs/bwa-0.7.12/bwa
BLAST = ${HOME}/programs/ncbi-blast-2.2.30+/bin/blastn
LAST = ${HOME}/programs/last-cea0cca54384

all: bam_idx bamBWA_idx

# Prog
bam_idx: results.bam
	samtools index $<

results.xml.gz: results.bam

results.bam: $(addprefix ${BINDIR}, blast2bam fastqParser) db.dict ${FASTQ} $(addsuffix .nin, ${DB})
	$(addprefix ${BINDIR}, fastqParser) ${FASTQ} | ${BLAST} -db ${DB} -num_threads 4 -word_size 8 -reward 1 -penalty -1 -gapopen 1 -gapextend 2 -outfmt 5 | \
	tee results.xml | $(addprefix ${BINDIR}, blast2bam) -W 40 -z -R '@RG	ID:foo	SM:sample' - db.dict ${FASTQ} | \
	samtools view -Sub - | samtools sort - $(basename $@) > $@
	gzip --best -f results.xml

$(addprefix ${BINDIR}, blast2bam fastqParser):
	(cd ../src; ${MAKE})

db.dict: ${DB}
	picard-tools CreateSequenceDictionary R=$< O=$@

$(addsuffix .nin,${DB}): ${DB}
	makeblastdb -in $< -dbtype 'nucl'

# BWA
bamBWA_idx: resultsBWA.bam
	samtools index $<

resultsBWA.bam: ${FASTQ} $(addsuffix .bwt, ${DB})
	${BWA} mem -t 4 -x 'ont2d' ${DB} ${FASTQ} | samtools view -Sub - | samtools sort - $(basename $@) > $@

$(addsuffix .bwt,${DB}): ${DB}
	${BWA} index ${DB}

# LAST
bamLAST_idx: resultsLASTSorted.bam
	samtools index $<

resultsLAST.bam: lastdb ${FASTQ} db.dict
	gunzip -c ${FASTQ} | ${LAST}/src/lastal -Q 1 $< | ${LAST}/scripts/maf-convert -f db.dict sam | \
	samtools view -Sub - | samtools sort - $(basename $@) > $@

lastdb: ${DB}
	${LAST}/src/lastdb $@ ${DB}

# Clean
clean:
	rm -f *.sam *.bam *.bai *.xml.gz db.dict lastdb*
	
destroy: clean
	(cd ../src; ${MAKE} clean)
