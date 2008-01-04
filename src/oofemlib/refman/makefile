# $Header: /home/cvs/bp/oofem/oofemlib/doc/refman/makefile,v 1.1 2001/06/12 15:18:55 bp Exp $

SHELL      = /bin/sh
BUILD_DIR_HTML  = html
BUILD_DIR_LATEX = latex
INST_DIR = oofemlibrefman
DOCXX      = doxygen
RCP = scp

all: htmldoc

htmldoc:
	$(DOCXX) doxygen.cfg

psdoc: htmldoc
	cd $(BUILD_DIR_LATEX); latex refman.tex; latex refman.tex; dvips -o refman.ps refman
	cp $(BUILD_DIR_LATEX)/refman.ps ..

www-install: htmldoc
ifdef INSTALLDIR
	mv $(BUILD_DIR_HTML) $(INST_DIR)
	$(RCP) -r $(INST_DIR) $(INSTALLDIR)/resources/doc
	mv $(INST_DIR) $(BUILD_DIR_HTML)
else
	echo "INSTALLDIR not defined"
endif

clean:
	rm -fr $(BUILD_DIR_HTML) $(BUILD_DIR_LATEX)

