SHELL = /bin/sh
BACKUP_FILES = *.tex *.fig *.eps makefile
RCP = scp

html: ps
	if [ x$(ROBOT) != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	latex2html -dir html -mkdir elementlibmanual

ps:
	if [ x$(ROBOT) != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	rm -f *.toc *.lot *.aux
	latex \\nonstopmode\\input elementlibmanual
	latex \\nonstopmode\\input elementlibmanual
	dvips -o 	elementlibmanual.ps elementlibmanual

pdf:
	if [ x$(ROBOT) != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	rm -f *.toc *.lot *.aux
	pdflatex elementlibmanual
	pdflatex elementlibmanual


www-install: html
ifdef INSTALLDIR
	$(RCP) -r html $(INSTALLDIR)/resources/doc/elementlibmanual
else
	echo "INSTALLDIR not defined"
endif

clean:
	rm -fr html
	rm -f *.aux *.dvi *.lof *.log *.toc