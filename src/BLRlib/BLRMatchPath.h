//
// Created by Hadi Quesneville on 2019-01-02.
//

#ifndef TE_FINDER_BLRMATHPATH_H
#define TE_FINDER_BLRMATHPATH_H

#include <iostream>
#include <stdlib.h>
#include <SDGString.h>
#include <utility>
#include <map>
#include <unistd.h>
#include "RangePairSet.h"
#include "BLRMatchAlign.h"
#include "BLRJoinParameter.h"
#include "BLRBioSeqDB.h"
#include "SDGBioSeqDB.h"

class BLRMatchPath {

public:

    // Key is numQ
    typedef std::map<long, std::list<RangePairSet> > MapPath;

    MapPath::iterator begin(){ return map_path.begin();};
    MapPath::iterator end(){ return map_path.end();};

    MapPath::const_iterator begin() const { return map_path.begin();};
    MapPath::const_iterator end() const { return map_path.end();};


private:

    MapPath map_path;
    std::map<std::string,long> name2numQ,name2numS;
    std::map<long,std::string> num2nameQ,num2nameS;
    unsigned countseqS,countseqQ;
    unsigned count_path;

public:
    BLRMatchPath(void){
        countseqS=0,countseqQ=0;count_path=0;
    };
    void insert(RangePairSet& range);

    std::list<RangePairSet>& operator[](long k){return map_path[k];};

    SDGString numQ2name(long num) { return num2nameQ[num];};
    SDGString numS2name(long num) { return num2nameS[num];};

    void setName2NumQ(std::map<std::string,long> name2NumQ){name2numQ=name2NumQ;};
    void setName2NumS(std::map<std::string,long> name2NumS){name2numS=name2NumS;};
    void setNum2NameQ(std::map<long,std::string> num2NameQ){num2nameQ=num2NameQ;};
    void setNum2NameS(std::map<long,std::string> num2NameS){num2nameS=num2NameS;};

    std::map<std::string,long> getName2NumQ(void){return name2numQ;};
    std::map<long, std::string> getNum2NameQ(void){return num2nameQ;};
    std::map<std::string,long> getName2NumS(void){return name2numS;};
    std::map<long, std::string> getNum2NameS(void){return num2nameS;};

    unsigned getNbQseq(void) { return num2nameQ.size();};
    unsigned getNbSseq(void) { return num2nameS.size();};
    unsigned getNbMatchesInMapPath(void){
        unsigned nbMatches = 0;
        for (MapPath::iterator m = map_path.begin(); m != map_path.end(); m++)
            for (std::list<RangePairSet>::iterator i = m->second.begin(); i != m->second.end(); i++)
                nbMatches += i->getNbRangePairs();
        return nbMatches;
    };
    unsigned getNbDistinctPaths(void){
        unsigned nbPaths = 0;
        for (MapPath::iterator m = map_path.begin(); m != map_path.end(); m++)
            nbPaths += m->second.size();
        return nbPaths;
    };


    void read(const BLRJoinParameter& param, std::istream& streamName, int verbose=0);
    void load(const BLRJoinParameter& param, const SDGString& filename, int verbose=0){
        std::ifstream input_align(filename);
        read(param,input_align, verbose);
    };
    void setFromRpsList(const BLRJoinParameter& p, const std::list<RangePairSet>& rps_list, int verbose=0);
    std::list<RangePairSet> getRpsListFromMapPath(void);
    void clear(void){map_path.clear();};


    void write(std::ostream &out);
    void write(const SDGString &filename){
        std::ostringstream pathStream;
        write(pathStream);
        std::ofstream pathFile(filename);
        pathFile << pathStream.str();
    };
    void writeAttribute(std::ostream &out);
    void writeBED(const SDGString &filename);
    void writeBED(std::ostream &out);
    void writeGFF3(const SDGString &filename);
    void writeGFF3(std::ostream &out);
    void writeMatch(const BLRJoinParameter& para, const SDGString &filename, int verbose);
    void writeSeq(const BLRJoinParameter& para, const SDGString &filename, int verbose);

};


#endif //TE_FINDER_BLRMATHPATH_H
