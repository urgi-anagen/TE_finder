#
# Makefile for SDGlib
#

include ../global_mak.def

OBJ = SDGString.o \
      SDGBioSeq.o \
      SDGMemBioSeq.o \
      SDGSubBioSeq.o \
      SDGFastaBioSeq.o \
      SDGBioSeqDB.o \
      SDGFastaIstream.o \
      SDGFastaOstream.o \
      SDGError.o \
      SDGBitStr.o \
      test.o

all:  libsdg.a 

libsdg.a: $(OBJ) 
	rm -f $@
	$(AR) rsc $@ $(OBJ)

install:
	@echo "Nothing to copy"

clean:
	@echo "clean the directory '"`pwd`"'"
	@find . -name '*.[oa]' -exec rm {} \;
	@find . -name 'core' -exec rm {} \;
	@rm -f *.so.*

test: libsdg.a  test.o
	$(CXX)  $(CXXFLAG)  $(INCLUDE) test.o -L. -lsdg -lstdc++ -lutil -o $@

test2: libsdg.a  test2.o
	$(CXX)  $(CXXFLAG)  $(INCLUDE) test.o -L. -lsdg -lstdc++ -o $@
