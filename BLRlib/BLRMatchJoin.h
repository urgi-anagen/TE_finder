//
// Created by Hadi Quesneville on 2019-01-02.
//

#ifndef TE_FINDER_BLRMATCHJOIN_H
#define TE_FINDER_BLRMATCHJOIN_H

#include <iostream>
#include <stdlib.h>
#include <SDGString.h>
#include <utility>
#include <unistd.h>
#include "BLRMatchAlign.h"
#include "BLRMatchPath.h"
#include "FragAlign.h"
#include "BLRJoinParameter.h"
#include "Graph.h"


class BLRMatchJoin {

    BLRJoinParameter para;
    unsigned merge_overlap;

    void add_clean_path_all_S(std::list<RangePairSet> &rp_list,
            std::list<RangePairSet>::iterator iter,
            BLRMatchPath &map_path);
    void add_clean_path_same_S(std::list<RangePairSet> &rp_list,
            std::list<RangePairSet>::iterator iter, BLRMatchPath &map_path);

    void add_split_path(std::list<RangePairSet> &rp_list, std::list<RangePairSet>::iterator iter, BLRMatchPath &map_path);

    public:

    BLRMatchJoin(BLRJoinParameter p): para(p), merge_overlap(200) {}


    void join(BLRMatchAlign &map_align, BLRMatchPath &map_path, int verbose=0);
    void noJoin(BLRMatchAlign &map_align, BLRMatchPath &map_path, int verbose=0 );
    void merge(BLRMatchPath &map_path, int verbose=0);
    void clean_conflicts(BLRMatchPath &map_path, int verbose=0);
    void split(BLRMatchPath &map_path, int verbose=0);

};


#endif //TE_FINDER_BLRMATCHJOIN_H
