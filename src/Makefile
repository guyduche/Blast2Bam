
.PHONY: all clean
CC = gcc
CFLAGS = -c -g -Wall `xml2-config --cflags`
LDFLAGS = -lz `xml2-config --libs`
VPATH = ../bin/
EXEC = prog
OBJS = parseXML.o blastSam.o shortRead.o utils.o main.o

all: ${OBJS} 
	${MAKE} ${EXEC}

%.o: %.c
	$(CC) ${CFLAGS} -o ../bin/$@ $<

parseXML.c: schema2c.xsl schema.xml
	xsltproc --output $@ --stringparam fileType c $^
 
parseXML.h: schema2c.xsl schema.xml
	xsltproc --output $@ --stringparam fileType h $^ 

prog: ${OBJS}
	$(CC) -o ../bin/$@ $^  ${LDFLAGS}

clean:
	rm -f ../bin/prog ../bin/*.o parseXML.c parseXML.h

main.o: blastSam.h
blastSam.o: blastSam.h parseXML.h utils.h shortRead.h
parseXML.o: parseXML.h utils.h
shortRead.o: shortRead.h utils.h
utils.o: utils.h
