SHELL = /bin/sh
BACKUP_FILES = *.tex *.fig *.eps makefile
RCP = scp

html: ps
	if [ x$(ROBOT) != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	latex2html -dir html -mkdir oofemInput

ps:
	if [ x$(ROBOT) != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	rm -f *.toc *.lot *.aux
	latex \\nonstopmode\\input oofemInput
	latex \\nonstopmode\\input oofemInput
	dvips -o 	oofemInput.ps oofemInput

pdf:
	if [ x$(ROBOT) != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	rm -f *.toc *.lot *.aux
	pdflatex oofemInput
	pdflatex oofemInput


www-install: html pdf
ifdef INSTALLDIR
	$(RCP) -r html $(INSTALLDIR)/resources/doc/oofemInput
	$(RCP) -r oofemInput.ps $(INSTALLDIR)/resources/doc/oofemInput
	$(RCP) -r oofemInput.pdf $(INSTALLDIR)/resources/doc/oofemInput
else
	echo "INSTALLDIR not defined"
endif

clean:
	rm -fr html
	rm -f *.aux *.dvi *.lof *.log *.toc
