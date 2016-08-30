#
# Makefile for BLRlib
#

include ../global_mak.def

OBJ =   BLRBioSeqDB.o \
	BLRJoinParameter.o \
	BLRBlasterParameter.o \
	BLRMatcherParameter.o \
	BLRGrouperParameter.o \
	BLRMatchList.o \
	BLRMatchMap.o \
	BLRRangePairSetListWriter.o\
	BLRMatchMapLoader.o
	

all: libBLR.a 

install:
	@echo "Nothing to copy"
	
clean:
	@echo "clean the directory '"`pwd`"'"
	@find . -name '*.[oa]' -exec rm {} \;
	@find . -name '*~' -exec rm {} \;
	@find . -name 'test.P' -exec rm {} \;
	@find . -name 'test' -exec rm {} \;

libBLR.a: $(OBJ) 
	rm -f $@
	$(AR) rsc $@ $(OBJ)


test: ../DGElib/libDGE.a libBLR.a $(OBJTEST) test.o 
	$(CXX) $(CXXFLAG) $(INCLUDE) $(LIB) $(OBJTEST) test.cpp -o $@ $(LINK) -lcppunit



