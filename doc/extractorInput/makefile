SHELL = /bin/sh
ROBOT = /home/bp/oofem-doc/www/robot.pl
BACKUP_FILES = *.tex makefile
RCP = scp

www-doc: ps-doc
	$(ROBOT) < latex2html-init.shtml > .latex2html-init
	latex2html extractorInput
	mv extractorInput.ps extractorInput

ps-doc:
	$(ROBOT) < include.src.tex > include.tex
	latex extractorInput
	latex extractorInput
	dvips -o extractorInput.ps extractorInput

www-install: www-doc
ifdef INSTALLDIR
	$(RCP) -r extractorInput $(INSTALLDIR)/resources/doc
else
	echo "INSTALLDIR not defined"
endif

clean:
	rm -fr extractorInput
	rm -f *.aux *.dvi *.lof *.log *.toc
