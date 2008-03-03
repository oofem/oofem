SHELL = /bin/sh

www-doc: 
	cd elementlibmanual; make www-doc
	cd oofemInput; make www-doc
	cd paper-actap; make www-doc
	cd extractorInput; make www-doc
	cd matlibmanual; make www-doc
	cd programmer; make www-doc
	cd refman; make www-doc

ps-doc:
	cd elementlibmanual; make ps-doc
	cd oofemInput; make ps-doc
	cd paper-actap; make ps-doc
	cd extractorInput; make ps-doc
	cd matlibmanual; make ps-doc
	cd programmer; make ps-doc


www-install:
ifdef INSTALLDIR
	export INSTALLDIR
	cd elementlibmanual; make www-install
	cd oofemInput; make www-install
	cd paper-actap; make www-install
	cd extractorInput; make www-install
	cd matlibmanual; make www-install
	cd programmer; make www-install
	cd refman; make www-install
else
	echo "INSTALLDIR not defined"
endif


clean:
	cd elementlibmanual; make clean
	cd oofemInput; make clean
	cd paper-actap; make clean
	cd extractorInput; make clean	
	cd matlibmanual; make clean
	cd programmer; make clean
	cd refman; make clean
