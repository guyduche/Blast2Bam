.PHONY: all clean
CC = gcc
CFLAGS = -c -O3 -Wall `xml2-config --cflags`
LDFLAGS = -lz `xml2-config --libs`
BINDIR = ../bin
EXEC = blast2bam
OBJS = $(addsuffix .o, parseXML blastSam printSam shortRead utils main)

all: ${BINDIR} ${BINDIR}/${EXEC}

${BINDIR}:
	mkdir ../bin

%.o: %.c
	$(CC) ${CFLAGS} -o $@ $<

parseXML.c: schema2c.xsl schema.xml
	xsltproc --output $@ --stringparam fileType c $^
 
parseXML.h: schema2c.xsl schema.xml
	xsltproc --output $@ --stringparam fileType h $^ 

${BINDIR}/blast2bam: ${OBJS}
	$(CC) -o $@ $^ ${LDFLAGS}

clean:
	rm -f ${BINDIR}/${EXEC} *.o parseXML.c parseXML.h

main.o: blastSam.h utils.h
blastSam.o: blastSam.h utils.h
printSam.o: blastSam.h utils.h
parseXML.o: parseXML.h utils.h
shortRead.o: shortRead.h utils.h
utils.o: utils.h
