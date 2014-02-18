SHELL = /bin/bash
BACKUP_FILES = *.tex *.fig *.eps makefile
RCP = scp

html: pdf
	if [ x"$(ROBOT)" != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	for i in elementlibmanual-figure*.pdf; do convert -density 300 $$i `echo $$i | sed -e 's/\.pdf/\.png/g'`; done
#mkdir html
	mkdir -p html
	htlatex elementlibmanual.tex "elementlibmanual.cfg,html,index=2,3,next" "" "-dhtml/" 
	rm elementlibmanual-figure*
	rm *.html

#latex2html -local_icons -dir html -mkdir elementlibmanual

ps: pdf
	pdftops elementlibmanual.pdf


pdf:
	if [ x"$(ROBOT)" != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	rm -f *.toc *.lot *.aux
	pdflatex -shell-escape elementlibmanual
	pdflatex -shell-escape elementlibmanual


www-install: html ps
ifdef INSTALLDIR
	$(RCP) -r html $(INSTALLDIR)/resources/doc/elementlibmanual
	$(RCP) -r elementlibmanual.ps $(INSTALLDIR)/resources/doc/elementlibmanual
	$(RCP) -r elementlibmanual.pdf $(INSTALLDIR)/resources/doc/elementlibmanual
else
	echo "INSTALLDIR not defined"
endif

clean:
	rm -fr html
	rm -f *.aux *.dvi *.lof *.log *.toc
