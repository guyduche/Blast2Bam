.PHONY: show clean
BASENAME=master_thesis

show: $(BASENAME).pdf 
	evince $< &

$(BASENAME).pdf: $(BASENAME).tex bibliography.bib
	latex $<
	makeglossaries $(BASENAME)
	makeglossaries $(BASENAME)
	bibtex $(BASENAME)
	latex $<
	latex $<
	latex $<
	dvips $(basename $<).dvi
	ps2pdf $(basename $<).ps

clean:
	rm -f $(addprefix $(BASENAME), .pdf .ps .dvi .toc .ist .glo .acn .acr .alg .lof .bbl .blg .glg .gls .out .log) *.aux
