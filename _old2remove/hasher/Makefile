#
#  Makefile for Blaster
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

$(TARGET): hasher_main.o Hasher.h Hasher.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG) -DVERSION=$(VERSION) $(INCLUDE) $(LIB) hasher_main.o Hasher.o -o $@ $(LINK)
