SHELL = /bin/sh
BACKUP_FILES = *.tex *.fig *.eps makefile
RCP = scp

www-doc: ps-doc
	if [ x$(ROBOT) != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	latex2html oofemInput
	mv oofemInput.ps oofemInput

ps-doc:
	if [ x$(ROBOT) != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	latex \\nonstopmode\\input oofemInput
	latex \\nonstopmode\\input oofemInput
	dvips -o 	oofemInput.ps oofemInput

www-install: www-doc
ifdef INSTALLDIR
	$(RCP) -r oofemInput $(INSTALLDIR)/resources/doc
else
	echo "INSTALLDIR not defined"
endif

clean:
	rm -fr oofemInput
	rm -f *.aux *.dvi *.lof *.log *.toc
