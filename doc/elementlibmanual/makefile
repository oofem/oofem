SHELL = /bin/sh
BACKUP_FILES = *.tex *.fig *.eps makefile
RCP = scp

www-doc: ps-doc
	$(ROBOT) < latex2html-init.shtml > .latex2html-init
	latex2html elementlibmanual
	mv elementlibmanual.ps elementlibmanual

ps-doc:
	$(ROBOT) < include.src.tex > include.tex
	latex elementlibmanual
	latex elementlibmanual
	dvips -o 	elementlibmanual.ps elementlibmanual

www-install: www-doc
ifdef INSTALLDIR
	$(RCP) -r elementlibmanual $(INSTALLDIR)/resources/doc
else
	echo "INSTALLDIR not defined"
endif

clean:
	rm -fr elementlibmanual
	rm -f *.aux *.dvi *.lof *.log *.toc