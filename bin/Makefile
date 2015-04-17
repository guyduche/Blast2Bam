
.PHONY: all clean
CC = gcc
CFLAGS = -c -g -Wall `xml2-config --cflags`
LDFLAGS = -lz `xml2-config --libs`
VPATH = ../src/
EXEC = prog parseDict
OBJS = parseXML.o cigar.o shortRead.o main.o

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

parseDict: parseDict.o
	$(CC) -o $@ $<

clean:
	rm -f prog parseDict *.o ../src/parseXML.c ../src/parseXML.h

parseDict.o: macro.h
main.o: parseXML.h cigar.h macro.h shortRead.h
cigar.o: cigar.h parseXML.h macro.h
parseXML.o: parseXML.h macro.h
shortRead.o: shortRead.h macro.h
