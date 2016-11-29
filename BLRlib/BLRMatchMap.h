/***
 *
 * BLRMatchMap.h
 *
 ***/

#ifndef BLRMATCHMAP_H
#define BLRMATCHMAP_H

#include <iostream>
#include <stdlib.h>
#include <SDGString.h>
#include <utility>
#include <map>
#include <unistd.h>
#include "RangePair.h"
#include "RangePairSet.h"
#include "RangeMap.h"
#include "BlastMatch.h"
#include "FragAlign.h"
#include "RangeSeq.h"
#include "BLRBioSeqDB.h"
#include "BLRMatcherParameter.h"
#include "Graph.h"

typedef std::list<RangePairSet> RpsList;

class BLRMatchMap
{
	friend class BLRMatchMapLoader;

 public:
  struct Key : std::pair<long,long>
  {
      Key(long i=-1,long j=-1)
      {
	first=i; second=j;
      };
  };

  // Key is numQ-numS
  typedef std::map<Key, std::list<RangePair> > MapAlign;
  typedef std::map<Key, std::list<RangePairSet> > MapPath;


 private:
  MapAlign map_align;
  MapPath map_path;
  BLRJoinParameter *para;
  std::map<std::string,long> name2numQ,name2numS;
  std::map<long,std::string> num2nameQ,num2nameS;
  BLRBioSeqDB subjectCut_db,queryCut_db;
  SDGBioSeqDB query_db,subject_db;
  bool same_db;
  RpsList rpsList;

  void add_clean(std::list<RangePair>& rp_list,
		      std::list<RangePair>::iterator iter);
  void insert_path(RangePairSet& range);
  void add_clean_path_same_S(std::list<RangePairSet>& rp_list,
		      std::list<RangePairSet>::iterator iter);
  // move to public for test and debug
  /*
  void add_clean_path_all_S(std::list<RangePairSet>& rp_list,
			    std::list<RangePairSet>::iterator iter,
			    int verbose=0);
  */
  static void insert_path_static(MapPath mapPath, RangePairSet& rangePair);
  void writePathForMergedS(std::ostream &out, unsigned path_id, std::string query_name, std::list<RangePair> path);

 public:

  BLRMatchMap(void)
  {
	  para=new BLRJoinParameter();
	  same_db=false;
  };

  BLRMatchMap(BLRJoinParameter *p) :
  para(p)
  {
	  same_db=false;
	  if(para->getQuery()!="<not set>")
		  query_db.load(para->getQuery());
	  if(para->getBank()!="<not set>")
		  subject_db.load(para->getBank());
  };


  BLRMatchMap(const BLRMatchMap& m)
  {
	  map_align=m.map_align;
	  map_path=m.map_path;
	  para=new BLRJoinParameter(*m.para);
	  rpsList=m.rpsList;
	  name2numQ=m.name2numQ;
	  name2numS=m.name2numS;
	  num2nameQ=m.num2nameQ;
	  num2nameS=m.num2nameS;
	  subjectCut_db=m.subjectCut_db;
	  queryCut_db=m.queryCut_db;
	  query_db=m.query_db;
	  subject_db=m.subject_db;
  	  same_db=m.same_db;
  }

  MapAlign::iterator begin(){ return map_align.begin();};
  MapAlign::iterator end(){ return map_align.end();};

  MapPath::iterator path_begin(){ return map_path.begin();};
  MapPath::iterator path_end(){ return map_path.end();};


  const SDGBioSeqDB& getRefQueryDB(void) {return query_db;};
  const SDGBioSeqDB& getRefSubjectDB(void) {return subject_db;};

  SDGString numQ2name(long num) { return num2nameQ[num];};
  SDGString numS2name(long num) { return num2nameS[num];};

  unsigned getNbQseq(void) { return num2nameQ.size();};
  unsigned getNbSseq(void);
  void readAlign(std::istringstream& streamName, int verbose=0);
  void load(int verbose=0){ load(para->getMatchFileName(),verbose);};
  void readPath(std::istringstream& streamName, int verbose=0);
  void loadPath(int verbose=0){ loadPath(para->getPath_filename(),verbose);};
  void load(SDGString filename, int verbose=0);
  void clear(void){map_align.clear();map_path.clear();};
  void mapPath(bool joining=true, bool clean_before=false, bool clean_after=false, bool merged=false, int verbose=0);
  void mapPathJoinOnlyForTest(bool joining=true, bool clean_before=false, bool clean_after=false, int verbose=0);
  //void mapPathWithThreads(bool joining=true, bool clean_before=false, bool clean_after=false, bool merged=false, int verbose=0);
  void clean_conflicts(void);
  void clean_path(bool same_S=true, int verbose=0);
  void split_path(void);
  //void clean_self(void);
  void select(bool subject=true, bool clean_before=false, bool clean_after=false);
  void selectQregex(SDGString regex);
  void writePath(const SDGString& filename, std::list<RangePairSet>& rps_list, int verbose=0);
  void writeRpsList(std::list<RangePairSet>& rps_list,std::ostream& out);
  void writeRpsListAttribute(std::list<RangePairSet>& rps_list, std::ostream& out);
  void writeBED(const SDGString& filename, const std::list<RangePairSet>& rps_list, const SDGString& color, int verbose);
  void writeBED(std::ostream& out, const std::list<RangePairSet>& rps_list, const SDGString& color, int verbose);
  void writeMatch(const SDGString& filename, int verbose=0);
  RangeMap writeMap(const SDGString& filename, int verbose=0);
  void writeSeq(const RangeMap& matchmap, const SDGString& filename, int verbose=0);
  void writeMapAlign(std::ostream& out);
  void contigOverlap(void);
  unsigned getNbMatchesInMapAlign(void);
  unsigned getNbMatchesInMapPath(void);
  unsigned getNbDistinctPaths(void);
  MapPath& getMapPath(void){return map_path;};
  
  std::list<RangePairSet> getRpsListFromMapPath(void);
  std::list<RangePairSet> copyRpsListFromMapPath(void);
  std::list<RangePairSet> getRpsList(void) // Build list of RangePairSet
	{ return rpsList;};
  void setRpsList(std::list<RangePairSet>& rpsListArg)
        { rpsList = rpsListArg;}; 
  void setMapPath(MapPath& mapPathArg){map_path = mapPathArg;};
  MapAlign getMapAlign(void){return map_align;};
  void insert(RangePair& range);
  void insert(RangePairSet& range);
  void insertRangePairSetIntoMapPath(RangePairSet rangePairSet, BLRMatchMap::MapPath mapPath);
  void add_clean_path_all_S(std::list<RangePairSet>& rp_list,
			    std::list<RangePairSet>::iterator iter,
			    int verbose=0);
  static bool isOverlapFound_in_add_split_path(std::list<RangePairSet>::iterator iter, MapPath mapPath, double idTolerance, unsigned lenFilter);
  void computeScoreWithLength(std::list<RangePairSet>& rpsList);
  void computeScoreWithLength();
  void add_split_path(std::list<RangePairSet>& rp_list, std::list<RangePairSet>::iterator iter);

  void merge(int verbose=0);
  std::list<RangePairSet> mergeOnCluster(std::list<RangePairSet> rpsList, Graph<unsigned long>& graph, int verbose=0);
  Graph <unsigned long> clusterizeOverlapingRps(const std::list<RangePairSet>& rpsList,int verbose=0);

  void loadPath(SDGString filename, int verbose=0);
  // for unit test purpose

  void setSameDb(bool sameDb){same_db = sameDb;};
  bool isSameDb(){return same_db;};

  void setName2NumQ(std::map<std::string,long> name2NumQ){name2numQ=name2NumQ;};
  void setName2NumS(std::map<std::string,long> name2NumS){name2numS=name2NumS;};
  void setNum2NameQ(std::map<long,std::string> num2NameQ){num2nameQ=num2NameQ;};
  void setNum2NameS(std::map<long,std::string> num2NameS){num2nameS=num2NameS;};

  std::map<std::string,long> getName2NumQ(void){return name2numQ;};
  std::map<long, std::string> getNum2NameQ(void){return num2nameQ;};
  std::map<std::string,long> getName2NumS(void){return name2numS;};
  std::map<long, std::string> getNum2NameS(void){return num2nameS;};

  BLRJoinParameter* getParameter(){return para;};
  
  void mapPathJoinAndComputeScoreWithLengthOnly(bool joining, bool clean_before, bool clean_after, int verbose);
  // TODO remove	
  void viewMapPath();
};

#endif
