SUBDIRS := $(wildcard */.)

all: $(SUBDIRS)
$(SUBDIRS):
	echo $(SUBDIRS)
	$(MAKE) -C $@
.PHONY: all $(SUBDIRS)
