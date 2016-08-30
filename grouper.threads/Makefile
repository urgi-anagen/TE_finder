#
# Makefile for grouper
#

include ../global_mak.def

//OBJ = BLRGroup.o BLRMemIdx.o
OBJ = BLRGroupThreads.o BLRGrouperThreads.o BLRMemIdxBin.o BLRBinSearch.o

TARGET = grouperThreads$(VERSION)

all: $(TARGET)

install: all
	@for i in $(TARGET); do cp $$i ../bin; done

clean:
	@echo "clean the directory '"`pwd`"'"
	@find . -name '*.[oa]' -exec rm {} \;
	@find . -name '*~' -exec rm {} \;
	@rm -f $(TARGET)

	
grouperThreads$(VERSION): grouperThreads.o $(OBJ) BLRGrouperThreads.h BLRGroupThreads.h
	$(CXX) $(CXXFLAG)  $(INCLUDE) $(LIB) BLRGrouperThreads.o BLRGroupThreads.o BLRBinSearch.o BLRMemIdxBin.o  grouperThreads.o -o $@ $(LINK)
	
