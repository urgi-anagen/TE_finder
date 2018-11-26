#
# Makefile for matcher
#

include ../global_mak.def

TARGET = matcher$(VERSION)

all: $(TARGET)

install: all
	@for i in $(TARGET); do cp $$i ../bin; done

clean:
	@echo "clean the directory '"`pwd`"'"
	@find . -name '*.[oa]' -exec rm {} \;
	@find . -name '*~' -exec rm {} \;
	@find . -name 'test.P' -exec rm {} \;
	@find . -name 'test' -exec rm {} \;
	@rm -f $(TARGET)

$(TARGET): matcher.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG) $(INCLUDE) $(INCLUDE_PVM) $(LIB) matcher.o -o $@ $(LINK)

