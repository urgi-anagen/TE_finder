#
#  Makefile for grouper
#

include ../global_mak.def

OBJ=BLRGroup.o BLRMemIdx.o

OBJ_TEST= BLRGroup.o\
	Test_BLRGroupUtils.o\
	Test_BLRGroup.o

TARGET= grouper$(VERSION)

all: $(TARGET)

install: all
	for i in $(TARGET); do cp $$i ../bin; done

clean:
	@echo "clean the directory"
	@find . -name '*.[oa]' -exec rm {} \;
	@find . -name '*~' -exec rm {} \;
	@rm -f $(TARGET)

grouper$(VERSION): grouper.o $(OBJ) BLRGroup.h
	$(CXX) $(CXXFLAG) $(INCLUDE) $(LIB) $(OBJ) grouper.o -o $@ $(LINK)

test: $(TARGET)
	$(CXX) $(CXXFLAG) $(INCLUDE) $(LIB) Test_F_grouper.cpp test.cpp -o $@ $(LINK) -lcppunit

test_u: $(TARGET) $(OBJ_TEST)
	$(CXX) $(CXXFLAG) $(INCLUDE) $(LIB) $(OBJ_TEST) test.cpp -o $@ $(LINK) -lcppunit
