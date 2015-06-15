
.PHONY: all clean
CC = gcc
CFLAGS = -c -O3 -Wall `xml2-config --cflags`
LDFLAGS = -lz `xml2-config --libs`
BINDIR = ../bin/
EXEC = blast2bam fastqParser
OBJS = $(addsuffix .o, parseXML blastSam shortRead utils main)

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

$(addprefix ${BINDIR}, fastqParser): fastqParser.o shortRead.o utils.o
	$(CC) -o $@ $^ ${LDFLAGS}

clean:
	rm -f $(addprefix ${BINDIR}, ${EXEC}) *.o parseXML.c parseXML.h

main.o: blastSam.h
blastSam.o: blastSam.h parseXML.h utils.h shortRead.h
parseXML.o: parseXML.h utils.h
shortRead.o: shortRead.h utils.h
utils.o: utils.h
fastqParser.o: shortRead.h utils.h
