# Blast2Sam

convert NCBI blastn + XML to SAM/BAM


#Compilation

```
make
```

# Usage

```bash
cat blast.xml | bam2bam [options] > out.sam
```

# Example

```bash
bin/fastqParser ${FASTQ} |\
	blastn -db ${DB} -num_threads 4 -word_size 8 -reward 1 -penalty -1 -gapopen 1 -gapextend 2 -outfmt 5 | \
	bin/blast2bam -W 40 -z -R '@RG	ID:foo	SM:sample' - db.dict ${FASTQ} | \
	samtools view -Sub - | samtools sort - $(basename $@) > $@
```

#Authors

* Aurelien Guy-Duche
* Pierre Lindenbaum


## Contribute

- Issue Tracker: http://github.com/lindenb/jvarkit/issues
- Source Code: http://github.com/lindenb/jvarkit


#License

The project is licensed under the MIT2 license.


