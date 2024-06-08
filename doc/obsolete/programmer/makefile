SHELL = /bin/sh
BACKUP_FILES = *.tex *.eps *.sty makefile
RCP = scp

html: ps
	if [ x"$(ROBOT)" != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	latex2html -style "/oofeml2h.css" -local_icons -dir html -mkdir programmer

ps:

	if [ x"$(ROBOT)" != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	rm -f *.toc *.lot *.aux
	latex \\nonstopmode\\input programmer
	latex \\nonstopmode\\input programmer
	dvips -o 	programmer.ps programmer

pdf:
	if [ x"$(ROBOT)" != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	rm -f *.toc *.lot *.aux
	pdflatex programmer
	pdflatex programmer

www-install: html pdf
ifdef INSTALLDIR
	$(RCP) -r html $(INSTALLDIR)/resources/doc/programmer
	$(RCP) -r programmer.ps $(INSTALLDIR)/resources/doc/programmer
	$(RCP) -r programmer.pdf $(INSTALLDIR)/resources/doc/programmer
else
	echo "INSTALLDIR not defined"
endif

clean:
	rm -fr programmer
	rm -f *.aux *.dvi *.lof *.log *.toc

