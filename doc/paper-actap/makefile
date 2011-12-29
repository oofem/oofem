SHELL = /bin/sh
BACKUP_FILES = *.tex *.eps *.sty makefile
RCP = scp

html: ps
	if [ x"$(ROBOT)" != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	latex2html -local_icons -dir html -mkdir paper-actap

ps:
	if [ x"$(ROBOT)" != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	rm -f *.toc *.lot *.aux
	latex \\nonstopmode\\input paper-actap
	latex \\nonstopmode\\input paper-actap
	dvips -o 	paper-actap.ps paper-actap

pdf:
	if [ x"$(ROBOT)" != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	rm -f *.toc *.lot *.aux
	pdflatex paper-actap
	pdflatex paper-actap

www-install: html pdf
ifdef INSTALLDIR
	$(RCP) -r html $(INSTALLDIR)/resources/doc/paper-actap
	$(RCP) -r paper-actap.ps $(INSTALLDIR)/resources/doc/paper-actap
	$(RCP) -r paper-actap.pdf $(INSTALLDIR)/resources/doc/paper-actap
else
	echo "INSTALLDIR not defined"
endif

clean:
	rm -fr paper-actap
	rm -f *.aux *.dvi *.lof *.log *.toc

