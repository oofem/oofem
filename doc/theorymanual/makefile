SHELL = /bin/sh
BACKUP_FILES = *.tex *.eps *.sty makefile
RCP = scp

html: ps
	latex2html -style "/oofeml2h.css" -local_icons -dir html -mkdir theorymanual

ps:

	rm -f *.toc *.lot *.aux
	latex \\nonstopmode\\input theorymanual
	latex \\nonstopmode\\input theorymanual
	dvips -o 	theorymanual.ps theorymanual

pdf:
	rm -f *.toc *.lot *.aux
	pdflatex theorymanual
	pdflatex theorymanual

www-install: html pdf
ifdef INSTALLDIR
	$(RCP) -r html $(INSTALLDIR)/resources/doc/theorymanual
	$(RCP) -r theorymanual.ps $(INSTALLDIR)/resources/doc/theorymanual
	$(RCP) -r theorymanual.pdf $(INSTALLDIR)/resources/doc/theorymanual
else
	echo "INSTALLDIR not defined"
endif

clean:
	rm -fr theorymanual
	rm -f *.aux *.dvi *.lof *.log *.toc

