#
# Makefile for grouper unitary tests
#

include ../../global_mak.def

BASE = ../..

OBJ = Test_Duster.o

OBJ_DEP = ../Duster.o 

TARGET = test_u 

all: $(TARGET)

clean:
	@echo "clean the directory '"`pwd`"'"
	@find . -name '*.[oa]' -exec rm {} \;
	@find . -name '*~' -exec rm {} \;
	@rm -f $(TARGET)

test_u: $(OBJ) $(OBJ_DEP) Test_Duster.h ../HashDNASeq.h
	$(CXX) $(CXXFLAG) $(INCLUDE) $(LIB) -lcppunit $(OBJ_DEP) $(OBJ) test_u.cpp -o $@ $(LINK)


