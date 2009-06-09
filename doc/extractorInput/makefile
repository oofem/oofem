SHELL = /bin/sh
BACKUP_FILES = *.tex makefile
RCP = scp

www-doc: ps-doc
	if [ x$(ROBOT) != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	latex2html extractorInput
	mv extractorInput.ps extractorInput

ps-doc:
	if [ x$(ROBOT) != "x" ]; then $(ROBOT) < include.src.tex > include.tex; fi
	latex \\nonstopmode\\input extractorInput
	latex \\nonstopmode\\input extractorInput
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
