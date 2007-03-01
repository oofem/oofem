SHELL = /bin/sh
ROBOT = /home/bp/oofem-doc/www/robot.pl
BACKUP_FILES = *.tex *.eps *.sty makefile
RCP = scp

www-doc: ps-doc
	$(ROBOT) < latex2html-init.shtml > .latex2html-init
	latex2html programmer
	mv programmer.ps programmer

ps-doc:
	$(ROBOT) < include.src.tex > include.tex
	latex programmer
	latex programmer
	dvips -o 	programmer.ps programmer

www-install: www-doc
ifdef INSTALLDIR
	$(RCP) -r programmer $(INSTALLDIR)/resources/doc
else
	echo "INSTALLDIR not defined"
endif

clean:
	rm -fr programmer
	rm -f *.aux *.dvi *.lof *.log *.toc

