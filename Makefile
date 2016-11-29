#
# Makefile for TE_finder
#

TARGETS = SDGlib DGElib BLRlib blaster grouper.bin grouper.threads matcher tools hasher duster

all: $(TARGETS)
	for i in $(TARGETS); do cd $$i; $(MAKE); cd ..; done

test: $(TARGETS)
	for i in $(TARGETS); do cd $$i; $(MAKE) test; cd ..; done

install: all
	@if test -d bin; then rm -r bin; fi; mkdir bin
	@for i in $(TARGETS); do cd $$i; $(MAKE) install ; cd ..; done

clean:
	@for i in $(TARGETS); do cd $$i; $(MAKE) clean; cd ..; done
	@find . -name '*~' -exec rm {} \;
	@find . -name '*.bak' -exec rm {} \;
	@find . -name '\#*\#' -exec rm {} \;
	@find . -name '*.o' -exec rm {} \;
	@find . -name '*.a' -exec rm {} \;
	@find . -name '*.P' -exec rm {} \;
	@if test -d bin; then rm -r bin; fi
