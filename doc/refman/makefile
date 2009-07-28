# $Header: /home/cvs/bp/oofem/doc/refman/makefile,v 1.8 2003/04/23 14:19:22 bp Exp $

SHELL      = /bin/sh
BUILD_DIR_HTML  = html
BUILD_DIR_LATEX = latex
INST_DIR = oofemrefman
DOCXX      = doxygen
RCP = scp
ROBOT = /home/bp/oofem-doc/www/robot2.pl

all: html 

html:
	if [ -x "$(ROBOT)" ] ; then $(ROBOT) < footer.shtml > footer.html; fi
	$(DOCXX) doxygen.cfg

ps: html
	cd $(BUILD_DIR_LATEX); latex refman.tex; latex refman.tex; dvips -o refman.ps refman
	cp $(BUILD_DIR_LATEX)/refman.ps ..

www-install: html
ifdef INSTALLDIR
	mv $(BUILD_DIR_HTML) $(INST_DIR)
	$(RCP) -r $(INST_DIR) $(INSTALLDIR)/resources/doc
	mv $(INST_DIR) $(BUILD_DIR_HTML)
else
	echo "INSTALLDIR not defined"
endif

clean:
	rm -fr $(BUILD_DIR_HTML) $(BUILD_DIR_LATEX)

