#
# Makefile for Test_BLRlib
#

include ../../global_mak.def
BASE = ../..

OBJ =   ../BLRBioSeqDB.o \
	../BLRJoinParameter.o \
	../BLRBlasterParameter.o \
	../BLRMatcherParameter.o \
	../BLRGrouperParameter.o \
	../BLRMatchList.o \
	../BLRMatchMap.o \
	../BLRMatchMapLoader.o
	

OBJTEST= Test_BLRMatchMapUtils.o \
	 Test_BLRMatchMapLoadUtils.o \
	 Test_BLRMatchMap.o
	 #Test_BLRMatchMapLoad.o \

TARGET = test_u

all: $(TARGET)
	
clean:
	@echo "clean the directory '"`pwd`"'"
	@find . -name '*.[oa]' -exec rm {} \;
	@find . -name '*~' -exec rm {} \;
	@find . -name 'test.P' -exec rm {} \;
	@find . -name 'test' -exec rm {} \;


test_u: ../../DGElib/libDGE.a ../libBLR.a $(OBJTEST) test_u.o 
	$(CXX) $(CXXFLAG) $(INCLUDE) $(LIB) $(OBJTEST) test_u.cpp -o $@ $(LINK) -lcppunit
	./test_u



