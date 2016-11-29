#
# Makefile for matcher tests
#

include ../../global_mak.def

OBJ = test.o

TARGET = test

test: $(OBJ) Test_F_matcher.h
	$(CXX) $(CXXFLAG) -I. -I../../SDGlib -I../../DGElib -L. -L../../DGElib -lcppunit $(OBJ) Test_F_matcher.cpp -o $@ -lDGE

clean:
	@echo "clean the directory '"`pwd`"'"
	@find . -name '*.[oa]' -exec rm {} \;
	@find . -name '*~' -exec rm {} \;
	@rm -f $(TARGET)