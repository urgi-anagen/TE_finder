#
# Makefile for tools
#

include ../global_mak.def

TARGETS = cutterDB \
			itrsearch \
			ltrsearch \
			TRsearch \
			polyAtail \
			TSDsearch \
			map2db \
			mapdist \
			mapsize \
			map2flank53 \
			map2flank5 \
			map2flank3 \
			mapview \
			mapOp \
			refalign \
			statalign \
			fastlalign \
			hrepeat \
			halign \
			hsearch \
			SWalign \
			NWalign \
			orienter \
			setnum2id \
			pathnum2id \
			align2piler \
			align2recon \
			filterAlign \
			tabnum2id \
			rpt_map

OBJTEST= Test_orienter.o

all: $(TARGETS)

install: all
	@for i in $(TARGETS); do cp $$i ../bin; done
	@cp rename.sh ../bin

clean:
	@echo "clean the directory '"`pwd`"'"
	@find . -name '*.[oa]' -exec rm {} \;
	@find . -name '*~' -exec rm {} \;

cutterDB:   cutterDB.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   cutterDB.o -o $@ $(LINK)

itrsearch:   itrsearch.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   itrsearch.o -o $@ $(LINK)

ltrsearch:   ltrsearch.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   ltrsearch.o -o $@ $(LINK)

TRsearch:   TRsearch.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   TRsearch.o -o $@ $(LINK)

TSDsearch:   TSDsearch.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   TSDsearch.o -o $@ $(LINK)
polyAtail:   polyAtail.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   polyAtail.o -o $@ $(LINK)

map2db:   map2db.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   map2db.o -o $@ $(LINK)

map2flank53:   map2flank53.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   map2flank53.o -o $@ $(LINK)
map2flank5:   map2flank5.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   map2flank5.o -o $@ $(LINK)
map2flank3:   map2flank3.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   map2flank3.o -o $@ $(LINK)

mapview:   mapview.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   mapview.o -o $@ $(LINK)

mapdist:   mapdist.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   mapdist.o -o $@ $(LINK)
	
mapsize:   mapsize.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   mapsize.o -o $@ $(LINK)

mapOp:   mapOp.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   mapOp.o -o $@ $(LINK)

refalign:   refalign.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   refalign.o -o $@ $(LINK)

statalign:   statalign.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   statalign.o -o $@ $(LINK)

orienter:   orienter.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   orienter.o -o $@ $(LINK)

orienter-debug:   orienter.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG) -g -o 0  $(INCLUDE) $(LIB)   orienter.o -o $@ $(LINK)

fastlalign:   fastlalign.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   fastlalign.o -o $@ $(LINK)

hrepeat:   hrepeat.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   hrepeat.o -o $@ $(LINK)

halign:   halign.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   halign.o -o $@ $(LINK)

hsearch:   hsearch.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   hsearch.o -o $@ $(LINK)

rpt_map:   rpt_map.o $(OBJ_LIB)
#	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   rpt_map.o -o $@ $(LINK)
	$(CXX)  $(INCLUDE) $(LIB)   rpt_map.o -o $@ $(LINK)

SWalign:   SWalign.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   SWalign.o -o $@ $(LINK)

NWalign:   NWalign.o $(OBJ_LIB)
	$(CXX) $(CXXFLAG)   $(INCLUDE) $(LIB)   NWalign.o -o $@ $(LINK)

setnum2id: setnum2id.o $(OBJ_LIB)
	$(CXX) $(CXXFLAGS)  $(INCLUDE) $(LIB) $@.cpp -o $@ $(LINK)

pathnum2id: pathnum2id.o $(OBJ_LIB)
	$(CXX) $(CXXFLAGS)  $(INCLUDE) $(LIB) $@.cpp -o $@ $(LINK)

align2piler: align2piler.o $(OBJ_LIB)
	$(CXX) $(CXXFLAGS)  $(INCLUDE) $(LIB) $@.cpp -o $@ $(LINK)

align2recon: align2recon.o $(OBJ_LIB)
	$(CXX) $(CXXFLAGS)  $(INCLUDE) $(LIB) $@.cpp -o $@ $(LINK)

filterAlign: filterAlign.o $(OBJ_LIB)
	$(CXX) $(CXXFLAGS)  $(INCLUDE) $(LIB) $@.cpp -o $@ $(LINK)

tabnum2id: tabnum2id.o $(OBJ_LIB)
	$(CXX) $(CXXFLAGS)  $(INCLUDE) $(LIB) $@.cpp -o $@ $(LINK)

test-orienter: orienter $(OBJTEST)
	$(CXX) $(CXXFLAG) $(INCLUDE) $(LIB) Test_orienter.cpp test.cpp -o $@ $(LINK) -lcppunit

