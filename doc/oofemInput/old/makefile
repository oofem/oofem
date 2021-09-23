SHELL = /bin/sh
BACKUP_FILES = *.tex *.fig *.eps makefile
RCP = scp

html: ps
	latex2html -style "/oofeml2h.css" -local_icons -dir html -mkdir oofemInput  -init_file "latex2html-init"

ps:
	rm -f *.toc *.lot *.aux
	latex \\nonstopmode\\input oofemInput
	latex \\nonstopmode\\input oofemInput
	dvips -o 	oofemInput.ps oofemInput

pdf:
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
