#
#  Makefile for Hasher
#

include ../global_mak.def


TARGET = hasher$(VERSION)

all: $(TARGET)

install: all
	@for i in $(TARGET); do cp $$i ../bin; done

clean:
	@echo "clean the directory '"`pwd`"'"
	@find . -name '*.[oa]' -exec rm {} \;
	@find . -name '*~' -exec rm {} \;
	@rm -f $(TARGET)

$(TARGET): hasher_main.o Hasher.h HashDNASeq.h Hasher.o ../duster/Duster.h $(OBJ_LIB)
	$(CXX) $(CXXFLAG) -DVERSION=$(VERSION) $(INCLUDE) $(LIB) hasher_main.o Hasher.o ../duster/Duster.o -o $@ $(LINK)
