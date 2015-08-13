.PHONY: clean destroy bam_idx bamBWA_idx
.SHELL = /bin/bash
BINDIR = ../../bin/
SRCDIR = ../../src/
FASTQ = test.fastq.gz
DB = ecoli.fasta

all: bam_idx bamBWA_idx

# Blast2Bam
bam_idx: results.bam
	samtools index $<

results.bam: $(addprefix ${BINDIR}, blast2bam fastq2fasta) db.dict ${FASTQ} $(addsuffix .nin, ${DB})
	$(addprefix ${BINDIR}, fastq2fasta) ${FASTQ} | \
	blastn -db ${DB} -num_threads 4 -word_size 8 -reward 1 -penalty -1 -gapopen 1 -gapextend 2 -outfmt 5 | \
	$(addprefix ${BINDIR}, blast2bam) -R '@RG\tID:foo\tSM:sample' - db.dict ${FASTQ} | \
	samtools view -Sub -F 0xF00 - | samtools sort - $(basename $@) > $@

$(addprefix ${BINDIR}, blast2bam fastq2fasta):
	(cd ${SRCDIR}; ${MAKE})

db.dict: ${DB}
	picard-tools CreateSequenceDictionary R=$< O=$@

$(addsuffix .nin, ${DB}): ${DB}
	makeblastdb -in $< -dbtype 'nucl'

# BWA
bamBWA_idx: resultsBWA.bam
	samtools index $<

resultsBWA.bam: ${FASTQ} $(addsuffix .bwt, ${DB})
	bwa mem -t 4 -x 'pacbio' ${DB} ${FASTQ} | samtools view -Sub -F 0xF00 - | samtools sort - $(basename $@) > $@

$(addsuffix .bwt, ${DB}): ${DB}
	bwa index ${DB}

# Clean
clean:
	rm -f *.bam *.bai db.dict *.amb *.ann *.bwt *.nhr *.nin *.nsq *.pac *.sa
	
destroy: clean
	(cd ${SRCDIR}; ${MAKE} clean)
