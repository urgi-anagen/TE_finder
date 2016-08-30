#
# Makefile for DGElib
#

include ../global_mak.def

OBJ = BlastMatch.o\
      Range.o\
      RangeAlign.o\
      RangeSeq.o\
      RangePair.o\
      RangePairSet.o\
      RangeMap.o\
      Galign.o\
      Lalign.o\
      FastLalign.o\
      ExtAlign.o\
      ItrAlign.o\
      FastExtAlign.o\
      HashRepeat.o\
      HashAlign.o\
      HashSearch.o\
      Cutter.o \
      Graph.o\
      Alea.o\
      FragAlign.o \
      FileUtils.o \
      HashAlignClone.o\
      DiagClone.o
#     FragAlignThreads.o \
#      RangeGroups.o\
#     ThreadPool.o
#      options.o
#      gnuplot_i.o	
#      RNG.o\
#      ACG.o\
#      Random.o\
#      SkipList.o\

TEST_OBJ = Test_FileUtils.o\
	Test_Range.o\
	Test_RangeAlignSet.o\
	Test_RangePair.o\
	Test_RangePairSet.o\
  	Test_RangePairSetUtils.o \
  	Test_HashAlign.o \
  	Test_HashAlignClone.o \
  	Test_DiagClone.o 

all: libDGE.a

install:
	@echo "Nothing to copy"

clean:
	@echo "clean the directory '"`pwd`"'"
	@find . -name '*.[oa]' -exec rm {} \;
	@find . -name '*~' -exec rm {} \;

libDGE.a: $(OBJ) 
	rm -f $@
	$(AR) rsc $@ $(OBJ)

test: libDGE.a $(TEST_OBJ) test.o
	$(CXX) $(CXXFLAG) $(INCLUDE) $(LIB) $(TEST_OBJ) test.cpp -o $@ $(LINK) -lcppunit
