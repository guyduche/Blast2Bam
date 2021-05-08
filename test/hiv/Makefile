.PHONY: clean destroy bam_idx bamBWA_idx
.SHELL = /bin/bash
BINDIR = ../../bin
SRCDIR = ../../src
FF = fastq2fasta
FASTQ = test_1.fastq.gz test_2.fastq.gz
DB = hiv.fasta

all: bam_idx bamBWA_idx

# Blast2Bam
bam_idx: results.bam
	samtools index $<

results.bam: ${BINDIR}/blast2bam ${FASTQ} ${DB}.nin
	${FF} ${FASTQ} | \
	blastn -db ${DB} -num_threads 4 -word_size 8 -outfmt 5 | \
	${BINDIR}/blast2bam -R '@RG\tID:foo\tSM:sample' - ${DB} ${FASTQ} | \
	samtools view -Sub -F 0xF00 - | samtools sort - > $@

${BINDIR}/blast2bam:
	(cd ${SRCDIR}; ${MAKE})

${DB}.nin: ${DB}
	makeblastdb -in $< -dbtype 'nucl'

# BWA
bamBWA_idx: resultsBWA.bam
	samtools index $<

resultsBWA.bam: ${FASTQ} ${DB}.bwt
	bwa mem -t 4 ${DB} ${FASTQ} | samtools view -Sub -F 0xF00 - | samtools sort - > $@

${DB}.bwt: ${DB}
	bwa index ${DB}

# Clean
clean:
	rm -f *.bam *.bai *.amb *.ann *.bwt *.nhr *.nin *.nsq *.pac *.sa
	
destroy: clean
	(cd ${SRCDIR}; ${MAKE} clean)
