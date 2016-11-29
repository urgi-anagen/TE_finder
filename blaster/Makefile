#
#  Makefile for blaster
#

include ../global_mak.def

TARGET = blaster$(VERSION)

all: $(TARGET)

install: all
	@for i in $(TARGET); do cp $$i ../bin; done

clean:
	@echo "clean the directory '"`pwd`"'"
	@find . -name '*.[oa]' -exec rm {} \;
	@find . -name '*~' -exec rm {} \;
	@rm -f $(TARGET)

$(TARGET): blaster2.o BLRNCBIBlast.o BLRNCBIBlastPlus.o BLRWUBlast.o BLRBlaster.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG) -DVERSION=$(VERSION) $(INCLUDE) $(LIB) BLRNCBIBlast.o BLRNCBIBlastPlus.o BLRWUBlast.o BLRBlaster.o blaster2.o -o $@ $(LINK)
