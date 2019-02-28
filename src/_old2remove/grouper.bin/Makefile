#
# Makefile for grouper
#

include ../global_mak.def

OBJ = BLRGroup.o BLRMemIdx.o

TARGET = grouper$(VERSION) 

all: $(TARGET)

install: all
	@for i in $(TARGET); do cp $$i ../bin; done

clean:
	@echo "clean the directory '"`pwd`"'"
	@find . -name '*.[oa]' -exec rm {} \;
	@find . -name '*~' -exec rm {} \;
	@rm -f $(TARGET)

grouper$(VERSION): grouper.o $(OBJ) BLRGroup.h
	$(CXX) $(CXXFLAG) $(INCLUDE) $(LIB) BLRGroup.o BLRMemIdx.o grouper.o -o $@ $(LINK)


