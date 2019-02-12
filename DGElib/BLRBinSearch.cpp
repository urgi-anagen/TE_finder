/**
 * \file BLRBinSearch.cpp
 */

#include "BLRBinSearch.h"

void threadedSearch(BLRBinSearch* bs, unsigned start, unsigned end, std::vector<unsigned>& vm)
{
   //std::mutex m;
   for( unsigned bin_lvl=bs->min_lvl; bin_lvl<=bs->max_lvl; bin_lvl++)
   {
        unsigned start_bin=bs->getBin(start,bin_lvl);
        unsigned end_bin=bs->getBin(end,bin_lvl);
        for(unsigned bin=start_bin; bin<=end_bin; bin++)
        {
            unsigned idx=unsigned((bin_lvl-bs->min_lvl+1)*bs->nb_bin_per_lvl)+bin;
            //std::unique_lock<std::mutex> lk(m);
            std::list<BLRBinSearch::coord>& lcoord=bs->mcoord[idx];
            //lk.unlock();
            for(std::list<BLRBinSearch::coord>::iterator i=lcoord.begin();i!=lcoord.end();i++)
                if((i->start>=start && i->start<=end)
                        ||(i->end>=start && i->end<=end)
                        ||(i->start<=start && i->end>=end))
                    vm.push_back(i->id);
        }
   }
};
