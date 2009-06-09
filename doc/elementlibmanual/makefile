SHELL = /bin/sh
BACKUP_FILES = *.tex *.fig *.eps makefile
RCP = scp

www-doc: ps-doc
	if [ x$(ROBOT) != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	latex2html elementlibmanual
	mv elementlibmanual.ps elementlibmanual

ps-doc:
	if [ x$(ROBOT) != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	latex \\nonstopmode\\input elementlibmanual
	latex \\nonstopmode\\input elementlibmanual
	dvips -o 	elementlibmanual.ps elementlibmanual

www-install: www-doc
ifdef INSTALLDIR
	$(RCP) -r elementlibmanual $(INSTALLDIR)/resources/doc
else
	echo "INSTALLDIR not defined"
endif

clean:
	rm -fr elementlibmanual
	rm -f *.aux *.dvi *.lof *.log *.toc