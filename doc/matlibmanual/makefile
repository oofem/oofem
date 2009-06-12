SHELL = /bin/sh
BACKUP_FILES = *.tex makefile
RCP = scp

www: ps
	if [ x$(ROBOT) != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	latex2html matlibmanual
	mv matlibmanual.ps matlibmanual

ps:
	if [ x$(ROBOT) != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	latex \\nonstopmode\\input matlibmanual
	latex \\nonstopmode\\input matlibmanual
	dvips -o 	matlibmanual.ps matlibmanual

pdf:
	if [ x$(ROBOT) != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	rm *.toc *.lot *.aux
	pdflatex matlibmanual


www-install: www-doc
ifdef INSTALLDIR
	$(RCP) -r matlibmanual $(INSTALLDIR)/resources/doc
else
	echo "INSTALLDIR not defined"
endif

clean:
	rm -fr matlibmanual
	rm -f *.aux *.dvi *.lof *.log *.toc