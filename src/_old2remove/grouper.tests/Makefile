#
# Makefile for grouper
#

include ../global_mak.def

OBJ = BLRGroup.o BLRMemIdx.o BLRGroupEBV.o BLRMemIdxEBV.o BLRMemIdxExhaustSearch.o BLRGroupExhaustSearch.o  \
BLRGroupBin.o BLRMemIdxBin.o BLRGroupVecBin.o BLRMemIdxVecBin.o BLRBinSearch.o

TARGET = grouper$(VERSION) grouperEBV$(VERSION) grouperBin$(VERSION) grouperVecBin$(VERSION) grouperExhaustSearch$(VERSION)

all: $(TARGET)

install: all
	@for i in $(TARGET); do cp $$i ../bin; done

clean:
	@echo "clean the directory '"`pwd`"'"
	@find . -name '*.[oa]' -exec rm {} \;
	@find . -name '*~' -exec rm {} \;
	@rm -f $(TARGET)

grouper$(VERSION): grouper.o $(OBJ) BLRGroup.h
	$(CXX) $(CXXFLAG) $(INCLUDE) $(LIB) BLRGroup.o BLRBinSearch.o BLRMemIdx.o grouper.o -o $@ $(LINK)


grouperEBV$(VERSION): grouperEBV.o $(OBJ) BLRGroupEBV.h
	$(CXX) $(CXXFLAG) $(INCLUDE) $(LIB) BLRGroupEBV.o BLRBinSearch.o BLRMemIdxEBV.o grouperEBV.o -o $@ $(LINK)


grouperBin$(VERSION): grouperBin.o $(OBJ) BLRGroupBin.h
	$(CXX) $(CXXFLAG)  $(INCLUDE) $(LIB) BLRGroupBin.o BLRBinSearch.o BLRMemIdxBin.o grouperBin.o -o $@ $(LINK)

grouperVecBin$(VERSION): grouperVecBin.o $(OBJ) BLRGroupVecBin.h
	$(CXX) $(CXXFLAG)  $(INCLUDE) $(LIB) BLRGroupVecBin.o BLRBinSearch.o BLRMemIdxVecBin.o grouperVecBin.o -o $@ $(LINK)
	
grouperExhaustSearch$(VERSION): grouperExhaustSearch.o $(OBJ) BLRGroupExhaustSearch.h
	$(CXX) $(CXXFLAG)  $(INCLUDE) $(LIB) BLRGroupExhaustSearch.o BLRBinSearch.o BLRMemIdxExhaustSearch.o grouperExhaustSearch.o -o $@ $(LINK)