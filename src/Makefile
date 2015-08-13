.PHONY: all clean
CC = gcc
CFLAGS = -c -O3 -Wall `xml2-config --cflags`
LDFLAGS = -lz `xml2-config --libs`
BINDIR = ../bin/
EXEC = blast2bam fastq2fasta
OBJS = $(addsuffix .o, parseXML blastSam printSam shortRead utils main)

all: ${BINDIR} $(addprefix ${BINDIR}, ${EXEC})

${BINDIR}:
	mkdir ../bin

%.o: %.c
	$(CC) ${CFLAGS} -o $@ $<

parseXML.c: schema2c.xsl schema.xml
	xsltproc --output $@ --stringparam fileType c $^
 
parseXML.h: schema2c.xsl schema.xml
	xsltproc --output $@ --stringparam fileType h $^ 

$(addprefix ${BINDIR}, blast2bam): ${OBJS}
	$(CC) -o $@ $^ ${LDFLAGS}

$(addprefix ${BINDIR}, fastq2fasta): fastq2fasta.o shortRead.o utils.o
	$(CC) -o $@ $^ ${LDFLAGS}

clean:
	rm -f $(addprefix ${BINDIR}, ${EXEC}) *.o parseXML.c parseXML.h

main.o: blastSam.h utils.h
blastSam.o: blastSam.h utils.h
printSam.o: blastSam.h utils.h
parseXML.o: parseXML.h utils.h
shortRead.o: shortRead.h utils.h
utils.o: utils.h
fastq2fasta.o: shortRead.h utils.h
