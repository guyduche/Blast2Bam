.PHONY: clean destroy bam_idx bamBWA_idx
.SHELL = /bin/bash
BINDIR = ../../bin/
SRCDIR = ../../src/
FASTQ = test_1.fastq.gz test_2.fastq.gz
DB = hiv2.fasta

all: bam_idx bamBWA_idx

# Blast2Bam
bam_idx: results.bam
	samtools index $<

results.bam: $(addprefix ${BINDIR}, blast2bam fastq2fasta) ${FASTQ} $(addsuffix .nin, ${DB})
	$(addprefix ${BINDIR}, fastq2fasta) ${FASTQ} | \
	blastn -db ${DB} -num_threads 4 -word_size 8 -outfmt 5 | \
	$(addprefix ${BINDIR}, blast2bam) -R '@RG\tID:foo\tSM:sample' - ${DB} ${FASTQ} | \
	samtools view -Sub -F 0xF00 - | samtools sort - $(basename $@) > $@

$(addprefix ${BINDIR}, blast2bam fastq2fasta):
	(cd ${SRCDIR}; ${MAKE})

$(addsuffix .nin, ${DB}): ${DB}
	makeblastdb -in $< -dbtype 'nucl'

# BWA
bamBWA_idx: resultsBWA.bam
	samtools index $<

resultsBWA.bam: ${FASTQ} $(addsuffix .bwt, ${DB})
	bwa mem -t 4 ${DB} ${FASTQ} | samtools view -Sub -F 0xF00 - | samtools sort - $(basename $@) > $@

$(addsuffix .bwt, ${DB}): ${DB}
	bwa index ${DB}

# Clean
clean:
	rm -f *.bam *.bai *.amb *.ann *.bwt *.nhr *.nin *.nsq *.pac *.sa
	
destroy: clean
	(cd ${SRCDIR}; ${MAKE} clean)
