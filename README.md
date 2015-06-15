# Blast2Bam

Map sequences with NCBI blastn and output the results in SAM/BAM format.


# Compilation

```
make
```

# Usage

```bash
	bin/fastqParser fastq1.gz [fastq2.gz] | \
	blastn [options] -db ref.fa -outfmt 5 | \
	bin/blast2bam [options] - ref.dict fastq1.gz [fastq2.gz] > out.sam
```

# Example

```bash
	bin/fastqParser ${FASTQ} |\
	blastn -db ${DB} -num_threads 4 -word_size 8 -reward 1 -penalty -1 -gapopen 1 -gapextend 2 -outfmt 5 | \
	bin/blast2bam -W 40 -z -R '@RG	ID:foo	SM:sample' - db.dict ${FASTQ} | \
	samtools view -Sub - | samtools sort - out > out.bam
```

# Authors

- Aurélien Guy-Duché
- Pierre Lindenbaum


# Contribute

- Issue Tracker: https://github.com/guyduche/prog/issues
- Source Code: https://github.com/guyduche/prog


# License

The project is licensed under the MIT2 license.


