#
# cible : dependance1 dependance2 ...
# (tab)azdopoakzdpokazdpokazdpkoazd
# (tab)apzdâzdkâzdkp
# $@ nom de la cible
# $< nom de la PREMIERE dependance
# $^ noms de TOUTES LES DEPENDANCES

.PHONY: clean
CC=gcc
CFLAGS=-c -g -Wall `xml2-config --cflags`
LDFLAGS=-lz `xml2-config --libs`

results: prog refDict
	gunzip -c ${HOME}/data/Pichon/180714_rand.fq.gz | head -n 40000 |paste - - - - | cut -f 1,2 | cut -c2- | sed 's/^/>/' | tr "\t" "\n" | blastn -db ~/data/Pichon/db/db.fasta -query - -outfmt 5 | ./$< - ${HOME}/data/Pichon/180714_rand.fq.gz >> results
	gedit results

refDict: parseDict dbDict
	./$< db.dict > results

parseDict: parseDict.o
	$(CC) -o $@ $<

parseDict.o: parseDict.c macro.h
	$(CC) ${CFLAGS} -o $@ $<

dbDict: ${HOME}/data/Pichon/db/db.fasta
	picard-tools CreateSequenceDictionary R=$< O=db.dict

prog: parseXML.o cigar.o main.o
	$(CC) -o $@ $^  ${LDFLAGS}
	
main.o: main.c parseXML.h cigar.h macro.h kseq.h
	$(CC) ${CFLAGS} -o $@ $<

cigar.o: cigar.c cigar.h parseXML.h macro.h
	$(CC) ${CFLAGS} -o $@ $<

parseXML.o: parseXML.c parseXML.h macro.h
	$(CC) ${CFLAGS} -o $@ $<

parseXML.c : schema2c.xsl schema.xml
	xsltproc --output $@ --stringparam fileType c $^
 
parseXML.h : schema2c.xsl schema.xml
	xsltproc --output $@ --stringparam fileType h $^ 

clean:
	rm -f results prog parseDict *.o parseXML.c parseXML.h db.dict

