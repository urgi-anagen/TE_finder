#
# Makefile for matcher tests
#

include ../../global_mak.def

BASE = ../..

OBJ = Test_F_matcher.o

OBJ_DEP = 

TARGET = test_u

all: $(TARGET)

clean:
	@echo "clean the directory '"`pwd`"'"
	@find . -name '*.[oa]' -exec rm {} \;
	@find . -name '*~' -exec rm {} \;
	@rm -f $(TARGET)

test_u: $(OBJ) $(OBJ_DEP) Test_F_matcher.h
	ln -sf ../matcher$(VERSION) matcher
	$(CXX) $(CXXFLAG) $(INCLUDE) $(LIB) -lcppunit $(OBJ_DEP) $(OBJ) test_u.cpp -o $@ $(LINK)
	./test_u

