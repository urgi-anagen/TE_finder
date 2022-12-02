//
// Created by Hadi Quesneville on 2019-01-22.
//

#ifndef TE_FINDER_BLRMATCHERTHREADS_H
#define TE_FINDER_BLRMATCHERTHREADS_H

#include "BLRMatchPath.h"
#include "BLRMatcherThreadsParameter.h"
#include "BLRMatchJoin.h"

class BLRMatcherThreads {

private:

    BLRMatcherThreadsParameter* matcher_parameter;
    BLRMatchPath map_path;

public:


    BLRMatcherThreads( BLRMatcherThreadsParameter* para):
            matcher_parameter(para)
    {};

    std::list<RangePairSet> getRpsListFromMapPath(void){
        return map_path.getRpsListFromMapPath();}
    void process(const std::list<RangePair>& rp_list);
};


#endif //TE_FINDER_BLRMATCHERTHREADS_H
