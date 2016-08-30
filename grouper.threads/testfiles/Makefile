#
# Makefile for grouper unitary tests
#

include ../../global_mak.def

BASE = ../..

OBJ = Test_BLRGroup.o Test_BLRGroupUtils.o

OBJ_DEP = ../BLRGroupThreads.o ../BLRMemIdxBin.o

TARGET = test_u test_f

all: $(TARGET)

clean:
	@echo "clean the directory '"`pwd`"'"
	@find . -name '*.[oa]' -exec rm {} \;
	@find . -name '*~' -exec rm {} \;
	@rm -f $(TARGET)

test_u: $(OBJ) $(OBJ_DEP) Test_BLRGroup.h 
	$(CXX) $(CXXFLAG) $(INCLUDE) $(LIB) -lcppunit $(OBJ_DEP) $(OBJ) test_u.cpp -o $@ $(LINK)

test_f: $(OBJ) $(OBJ_DEP)
	$(CXX) $(CXXFLAG) $(INCLUDE) $(LIB) -lcppunit $(OBJ_DEP) $(OBJ) test_f.cpp -o $@ $(LINK)

