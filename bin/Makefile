
.PHONY: all clean
CC = gcc
CFLAGS = -c -g -Wall `xml2-config --cflags`
LDFLAGS = -lz `xml2-config --libs`
VPATH = ../src/
EXEC = prog
OBJS = parseXML.o blastSam.o shortRead.o main.o

all: parseXML.c parseXML.h 
	${MAKE} ${EXEC}

%.o: %.c
	$(CC) ${CFLAGS} -o $@ $<

parseXML.c: schema2c.xsl schema.xml
	xsltproc --output ../src/$@ --stringparam fileType c $^
 
parseXML.h: schema2c.xsl schema.xml
	xsltproc --output ../src/$@ --stringparam fileType h $^ 

prog: ${OBJS}
	$(CC) -o $@ $^  ${LDFLAGS}

clean:
	rm -f prog parseDict *.o ../src/parseXML.c ../src/parseXML.h

main.o: parseXML.h sam.h utils.h shortRead.h
blastSam.o: blastSam.h parseXML.h utils.h shortRead.h
parseXML.o: parseXML.h utils.h
shortRead.o: shortRead.h utils.h
