SHELL = /bin/sh
BACKUP_FILES = *.tex makefile
RCP = scp

html: ps
	if [ x"$(ROBOT)" != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	latex2html -style "/oofeml2h.css" -local_icons -dir html -mkdir extractorInput

ps:
	if [ x"$(ROBOT)" != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	rm -f *.toc *.lot *.aux
	latex \\nonstopmode\\input extractorInput
	latex \\nonstopmode\\input extractorInput
	dvips -o extractorInput.ps extractorInput

pdf:
	if [ x"$(ROBOT)" != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	rm -f *.toc *.lot *.aux
	pdflatex extractorInput
	pdflatex extractorInput


www-install: html pdf
ifdef INSTALLDIR
	$(RCP) -r html $(INSTALLDIR)/resources/doc/extractorInput
	$(RCP) -r extractorInput.ps $(INSTALLDIR)/resources/doc/extractorInput
	$(RCP) -r extractorInput.pdf $(INSTALLDIR)/resources/doc/extractorInput
else
	echo "INSTALLDIR not defined"
endif

clean:
	rm -fr extractorInput
	rm -f *.aux *.dvi *.lof *.log *.toc
