SHELL = /bin/sh
BACKUP_FILES = *.tex makefile
RCP = scp

html: ps
	latex2html -style "/oofeml2h.css" -local_icons -dir html -mkdir matlibmanual -init_file "latex2html-init"

ps:
	rm -f *.toc *.lot *.aux
	latex \\nonstopmode\\input matlibmanual
	latex \\nonstopmode\\input matlibmanual
	dvips -o 	matlibmanual.ps matlibmanual

pdf:
	rm -f *.toc *.lot *.aux
	pdflatex matlibmanual
	pdflatex matlibmanual


www-install: html pdf
ifdef INSTALLDIR
	$(RCP) -r html $(INSTALLDIR)/resources/doc/matlibmanual
	$(RCP) -r matlibmanual.ps $(INSTALLDIR)/resources/doc/matlibmanual
	$(RCP) -r matlibmanual.pdf $(INSTALLDIR)/resources/doc/matlibmanual
else
	echo "INSTALLDIR not defined"
endif

clean:
	rm -fr html
	rm -f *.aux *.dvi *.lof *.log *.toc
