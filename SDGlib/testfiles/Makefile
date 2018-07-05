#
# Makefile for grouper unitary tests
#

include ../../global_mak.def

BASE = ../..

OBJ = Test_SDGlib.o

OBJ_DEP = ../SDGBioSeq.o \
../SDGBioSeqDB.o \
../SDGBitStr.o \
../SDGError.o \
../SDGFastaBioSeq.o \
../SDGFastaIstream.o \
../SDGFastaOstream.o \
../SDGMemBioSeq.o \
../SDGString.o \
../SDGSubBioSeq.o

HEADERS_DEP = ../SDGBioSeq.h \
../SDGBioSeqDB.h \
../SDGBitStr.h \
../SDGError.h \
../SDGFastaBioSeq.h \
../SDGFastaIstream.h \
../SDGFastaOstream.h \
../SDGMemBioSeq.h \
../SDGString.h \
../SDGSubBioSeq.h


TARGET = test_u 

all: $(TARGET)

clean:
	@echo "clean the directory '"`pwd`"'"
	@find . -name '*.[oa]' -exec rm {} \;
	@find . -name '*~' -exec rm {} \;
	@rm -f $(TARGET)

	
test_u: $(OBJ) $(OBJ_DEP) $(HEARDERS_DEP)
	$(CXX) $(CXXFLAG) $(INCLUDE) $(LIB) -lcppunit $(OBJ_DEP) $(OBJ) test_u.cpp -o $@  $(LINK)


