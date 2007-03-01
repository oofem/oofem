SHELL = /bin/sh
ROBOT = /home/bp/oofem-doc/www/robot.pl
BACKUP_FILES = *.tex *.fig *.eps makefile
RCP = scp

www-doc: ps-doc
	$(ROBOT) < latex2html-init.shtml > .latex2html-init
	latex2html oofemInput
	mv oofemInput.ps oofemInput

ps-doc:
	$(ROBOT) < include.src.tex > include.tex
	latex oofemInput
	latex oofemInput
	dvips -o 	oofemInput.ps oofemInput

www-install: www-doc
ifdef INSTALLDIR
	$(RCP) -r oofemInput $(INSTALLDIR)/resources/doc
else
	echo "INSTALLDIR not defined"
endif

clean:
	rm -fr oofemInput
	rm -f *.aux *.dvi *.lof *.log *.toc
