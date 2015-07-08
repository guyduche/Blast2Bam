.PHONY: clean destroy
.SHELL = /bin/bash
DB = hiv.fasta
DB2 = hiv2.fasta

# hiv2.fasta contains hiv.fasta + a copy of hiv.fasta with a different name

# Tested with:
# BWA version 0.7.12
# SRA toolkit version 2.5.2-ubuntu64

testResultsExemple.txt: $(addprefix resultsBWA.test, $(addsuffix .sam, 1 2 3 4 5 6))
	$(foreach SAM, $^, \
		echo $(SAM) >> $@ ; \
		grep -Fw 'ERR656485.14' $(SAM) >> $@ ; \
		echo "\n" >> $@ ;)

# Test 1: single ref
resultsBWA.test1.sam: ERR656485_1.fastq.gz ERR656485_2.fastq.gz $(addsuffix .bwt, ${DB})
	bwa mem -t 4 ${DB} ERR656485_1.fastq.gz ERR656485_2.fastq.gz > $@

# Test 2: double ref
resultsBWA.test2.sam: ERR656485_1.fastq.gz ERR656485_2.fastq.gz $(addsuffix .bwt, ${DB2})
	bwa mem -t 4 ${DB2} ERR656485_1.fastq.gz ERR656485_2.fastq.gz > $@

# Test 3: single ref + interleaved mode
resultsBWA.test3.sam: ERR656485_1.test.fastq ERR656485_2.test.fastq $(addsuffix .bwt, ${DB})
	paste ERR656485_1.test.fastq ERR656485_2.test.fastq | tr "\t" "\n" | \
	bwa mem -p -t 4 ${DB} - > $@

# Test 4: double ref + interleaved mode
resultsBWA.test4.sam: ERR656485_1.test.fastq ERR656485_2.test.fastq $(addsuffix .bwt, ${DB2})
	paste ERR656485_1.test.fastq ERR656485_2.test.fastq | tr "\t" "\n" | \
	bwa mem -p -t 4 ${DB2} - > $@

# Test 5: single ref + interleaved mode + fastq shuffled
resultsBWA.test5.sam: ERR656485_1.test.fastq ERR656485_2.test.fastq $(addsuffix .bwt, ${DB})
	paste ERR656485_1.test.fastq ERR656485_2.test.fastq | \
	shuf | tr "\t" "\n" | \
	bwa mem -p -t 4 ${DB} - > $@

# Test 6: double ref + interleaved mode + fastq shuffled
resultsBWA.test6.sam: ERR656485_1.test.fastq ERR656485_2.test.fastq $(addsuffix .bwt, ${DB2})
	paste ERR656485_1.test.fastq ERR656485_2.test.fastq | \
	shuf | tr "\t" "\n" | \
	bwa mem -p -t 4 ${DB2} - > $@


# Temp files for tests
ERR656485_1.fastq.gz ERR656485_2.fastq.gz: ERR656485.sra
	fastq-dump --gzip --split-files $<

ERR656485.sra:
	curl "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/ERR656/ERR656485/ERR656485.sra" -o $@

$(addsuffix .bwt, ${DB}): ${DB}
	bwa index ${DB}

$(addsuffix .bwt, ${DB2}): ${DB2}
	bwa index ${DB2}
	
ERR656485_1.test.fastq: ERR656485_1.fastq.gz
	gunzip -c $< | paste - - - - > $@

ERR656485_2.test.fastq: ERR656485_2.fastq.gz
	gunzip -c $< | paste - - - - > $@

# Clean
clean:
	rm -f *.sam *.amb *.ann *.bwt *.pac *.sa *.test.fastq *.txt

destroy: clean
	rm -f *.gz *.sra
