SHELL = /bin/sh
BACKUP_FILES = *.tex *.eps *.sty makefile
RCP = scp

www-doc: ps-doc
	if [ x$(ROBOT) != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	latex2html paper-actap
	mv paper-actap.ps paper-actap

ps-doc:
	if [ x$(ROBOT) != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	latex \\nonstopmode\\input paper-actap
	latex \\nonstopmode\\input paper-actap
	dvips -o 	paper-actap.ps paper-actap

www-install: www-doc
ifdef INSTALLDIR
	$(RCP) -r paper-actap $(INSTALLDIR)/resources/doc
else
	echo "INSTALLDIR not defined"
endif

clean:
	rm -fr paper-actap
	rm -f *.aux *.dvi *.lof *.log *.toc

