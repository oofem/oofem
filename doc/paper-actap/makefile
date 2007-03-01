SHELL = /bin/sh
ROBOT = /home/bp/oofem-doc/www/robot.pl
BACKUP_FILES = *.tex *.eps *.sty makefile
RCP = scp

www-doc: ps-doc
	$(ROBOT) < latex2html-init.shtml > .latex2html-init
	latex2html paper-actap
	mv paper-actap.ps paper-actap

ps-doc:
	$(ROBOT) < include.src.tex > include.tex
	latex paper-actap
	latex paper-actap
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

