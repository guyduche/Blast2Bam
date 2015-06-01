.PHONY: clean destroy diff bam_idx bamBWA_idx
.SHELL=/bin/bash
VPATH = ../bin/
FASTQ = ${HOME}/data/nanopore/phage/ERR732557.fastq.gz
DB = ${HOME}/data/nanopore/phage/db/m13mp18.fasta
SRA = ${HOME}/data/HIV/ERR656485.sra

all: bam_idx bamBWA_idx

diff: results.sam resultsBWA.sam
	$(foreach S,$^,cut -f 1,3,4,6 $S | grep -vF '@SQ' | grep -F 'M13' | sed 's/	/~/' | LC_ALL=C sort -t '	' -k 1,1 -k 2,2  > $(addsuffix .tmp,$S); )
	join -t '	' -1 1 -2 1 $(addsuffix .tmp,$^) > same.read-chrom.txt
	join -t '	' -v 1 -1 1 -2 1 $(addsuffix .tmp,$^) > only.prog.txt
	join -t '	' -v 2 -1 1 -2 1 $(addsuffix .tmp,$^) > only.bwa.txt

# Prog
bam_idx: resultsSorted.bam
	samtools index $<

resultsSorted.bam: results.sam
	samtools view -Sb $< | samtools sort - $(basename $@)

results.xml.gz : resultsSorted.sam

results.sam : ./../bin/prog ./../bin/fastqParser db.dict ${FASTQ} $(addsuffix .nin,${DB})
	./../bin/fastqParser ${FASTQ} | blastn -db ${DB} -num_threads 4 -word_size 20 -outfmt 5 | tee results.xml |  ./../bin/prog - db.dict ${FASTQ} > $@
	gzip --best -f results.xml

./../bin/prog ./../bin/fastqParser:
	(cd ../src; ${MAKE} )

db.dict: ${DB}
	picard-tools CreateSequenceDictionary R=$< O=$@

$(addsuffix .nin,${DB}) : ${DB}
	makeblastdb -in $< -dbtype nucl

# BWA
bamBWA_idx: resultsBWASorted.bam
	samtools index $<

resultsBWASorted.bam: resultsBWA.sam
	samtools view -Sb $< | samtools sort - $(basename $@)

resultsBWA.sam: ${FASTQ} ${DB}
	bwa mem -k 20 ${DB} ${FASTQ} > $@

indexBWA: ${DB}
	bwa index ${DB}

# Fastq
fastq: ${SRA}
	(cd ${HOME}/data/HIV; fastq-dump --split-files --gzip ${SRA})

# Clean
clean:
	rm -f *.sam *.bam *.bai *.tmp *.txt *.xml.gz db.dict
	
destroy: clean
	(cd ../src; ${MAKE} clean)
