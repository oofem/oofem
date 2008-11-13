SHELL = /bin/sh
BACKUP_FILES = *.tex makefile
RCP = scp

www-doc: ps-doc
	$(ROBOT) < latex2html-init.shtml > .latex2html-init
	latex2html matlibmanual
	mv matlibmanual.ps matlibmanual

ps-doc:
	$(ROBOT) < include.src.tex > include.tex
	latex matlibmanual
	latex matlibmanual
	dvips -o 	matlibmanual.ps matlibmanual

www-install: www-doc
ifdef INSTALLDIR
	$(RCP) -r matlibmanual $(INSTALLDIR)/resources/doc
else
	echo "INSTALLDIR not defined"
endif

clean:
	rm -fr matlibmanual
	rm -f *.aux *.dvi *.lof *.log *.toc