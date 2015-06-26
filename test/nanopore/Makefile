.PHONY: clean destroy bam_idx bamBWA_idx
.SHELL=/bin/bash
BINDIR = ../../bin/
SRCDIR = ../../src/
FASTQ = ERR732557.fastq.gz
DB = m13mp18.fasta
BWA = ${HOME}/programs/bwa-0.7.12/bwa
BLAST = ${HOME}/programs/ncbi-blast-2.2.31+/bin

all: bam_idx bamBWA_idx

# Blast2Bam
bam_idx: results.bam
	samtools index $<

results.bam: $(addprefix ${BINDIR}, blast2bam fastq2fasta) db.dict ${FASTQ} $(addsuffix .nin, ${DB})
	$(addprefix ${BINDIR}, fastq2fasta) ${FASTQ} | \
	${BLAST}/blastn -db ${DB} -num_threads 4 -word_size 8 -reward 1 -penalty -1 -gapopen 1 -gapextend 2 -outfmt 5 | \
	$(addprefix ${BINDIR}, blast2bam) -W 40 -R '@RG	ID:foo	SM:sample' - db.dict ${FASTQ} | \
	samtools view -Sub - | samtools sort - $(basename $@) > $@

$(addprefix ${BINDIR}, blast2bam fastq2fasta):
	(cd ${SRCDIR}; ${MAKE})

db.dict: ${DB}
	picard-tools CreateSequenceDictionary R=$< O=$@

$(addsuffix .nin,${DB}): ${DB}
	${BLAST}/makeblastdb -in $< -dbtype 'nucl'

# BWA
bamBWA_idx: resultsBWA.bam
	samtools index $<

resultsBWA.bam: ${FASTQ} $(addsuffix .bwt, ${DB})
	${BWA} mem -t 4 -x 'ont2d' ${DB} ${FASTQ} | samtools view -Sub - | samtools sort - $(basename $@) > $@

$(addsuffix .bwt,${DB}): ${DB}
	${BWA} index ${DB}

# Clean
clean:
	rm -f *.sam *.bam *.bai db.dict *.amb *.ann *.bwt *.nhr *.nin *.nsq *.pac *.sa
	
destroy: clean
	(cd ${SRCDIR}; ${MAKE} clean)
