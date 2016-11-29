#
#  Makefile for Blaster
#

include ../global_mak.def


TARGET = duster$(VERSION)

all: $(TARGET)

install: all
	@for i in $(TARGET); do cp $$i ../bin; done

clean:
	@echo "clean the directory '"`pwd`"'"
	@find . -name '*.[oa]' -exec rm {} \;
	@find . -name '*~' -exec rm {} \;
	@rm -f $(TARGET)

$(TARGET): duster_main.o Duster.h HashDNASeq.h Duster.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG) -DVERSION=$(VERSION) $(INCLUDE) $(LIB) duster_main.o Duster.o -o $@ $(LINK)
