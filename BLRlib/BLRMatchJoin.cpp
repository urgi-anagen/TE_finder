//
// Created by Hadi Quesneville on 2019-01-02.
//

#include "BLRMatchJoin.h"

//------------------------------------------------------------------------------------------------------------
void BLRMatchJoin::join(BLRMatchAlign &map_align, BLRMatchPath &map_path, int verbose ) {
    if (verbose > 0)
        std::cout << "Join parameters: dist_pen=" << para.getDist_pen()
                       << " gap_pen=" << para.getGap_pen()
                       << " overlap=" << para.getOverlap() << std::__1::endl << std::__1::flush;

    unsigned count = 0;
    for (BLRMatchAlign::MapAlign::iterator m = map_align.begin(); m != map_align.end(); m++) {
        std::list<RangePairSet> &list = map_path[m->first.first];
        FragAlign fragAlign(para.getDist_pen(), 0, para.getGap_pen(),
                            para.getOverlap());
        list.splice(list.end(),fragAlign.join(m->second));
        for(std::list<RangePairSet>::iterator it=list.begin(); it!=list.end();it++)
            it->setId(count++);
        m->second.clear();
    }
    map_path.setName2NumQ(map_align.getName2NumQ());
    map_path.setName2NumS(map_align.getName2NumS());
    map_path.setNum2NameQ(map_align.getNum2NameQ());
    map_path.setNum2NameS(map_align.getNum2NameS());

    if (verbose > 0) {
        std::cout << "After join:\nnb of matches: " << map_path.getNbMatchesInMapPath() << std::__1::endl;
        std::cout << "nb of paths: " << map_path.getNbDistinctPaths() << std::__1::endl;
    }
}
//------------------------------------------------------------------------------------------------------------
void BLRMatchJoin::noJoin(BLRMatchAlign &map_align, BLRMatchPath &map_path, int verbose ) {
    unsigned count = 0;
    for (BLRMatchAlign::MapAlign::iterator m = map_align.begin(); m != map_align.end(); m++) {
        std::list<RangePairSet> &list = map_path[m->first.first];
        for (std::list<RangePair>::iterator i = m->second.begin();
             i != m->second.end(); i++) {
            RangePairSet rps(*i);
            rps.setId(count++);
            list.push_back(rps);
        }
        m->second.clear();
    }
    map_path.setName2NumQ(map_align.getName2NumQ());
    map_path.setName2NumS(map_align.getName2NumS());
    map_path.setNum2NameQ(map_align.getNum2NameQ());
    map_path.setNum2NameS(map_align.getNum2NameS());

    if (verbose > 0)
        std::cout << count << " matches" << std::endl;
}
//----------------------------------------------------------------------------
void BLRMatchJoin::add_clean_path_all_S(std::list<RangePairSet> &rp_list,
                                       std::list<RangePairSet>::iterator iter,
                                       BLRMatchPath &map_path)
//add a rangePairSet and post process it removing conflicting subjects
{
    bool found_over = false;
    bool atLeastOneOverlap = false;
    unsigned nbseqS = map_path.getNbSseq();
    //unsigned s=1;
    // start subject at -1 to consider merged data
    long s = -1;
    // for current query iter on each subject in map_path.
    while (s <= nbseqS) {
        // list of subject
        std::list<RangePairSet> &list
                = map_path[iter->getRangeQ().getNumChr()];
        for (std::list<RangePairSet>::iterator iter_list = list.begin();
             iter_list != list.end(); iter_list++) {
            if (RangePair::greaterScore(*iter_list, *iter)
                && iter->overlapQ(*iter_list))
                if (iter->diffQ(*iter_list))
                    found_over = true;
        }
        if (found_over) {
            atLeastOneOverlap = true;
            if (!iter->empty()
                && iter->getRangeQ().getLength() > para.getLenFilter()) {

                std::list<RangePairSet>::iterator it
                        = std::lower_bound(iter, rp_list.end(),
                                           *iter,
                                           RangePair::greaterScore);
                if (it == iter)
                    it++;
                rp_list.insert(it, *iter);


            }
        }

        if (!found_over) s++; else found_over = false;
    }
    if (!atLeastOneOverlap && iter->getRangeQ().getLength() > para.getLenFilter()
        && iter->getScore() > 0) {
        map_path.insert(*iter);
    }
}
//----------------------------------------------------------------------------
void BLRMatchJoin::clean_conflicts(BLRMatchPath &map_path, int verbose)
// removing conflicting subjects
{
    BLRMatchPath map_path_new;
    map_path_new.setName2NumQ(map_path.getName2NumQ());
    map_path_new.setName2NumS(map_path.getName2NumS());
    map_path_new.setNum2NameQ(map_path.getNum2NameQ());
    map_path_new.setNum2NameS(map_path.getNum2NameS());

    for (BLRMatchPath::MapPath::iterator m = map_path.begin(); m != map_path.end(); m++) {
        m->second.sort(RangePair::greaterScore);
        for (std::list<RangePairSet>::iterator i = m->second.begin();
             i != m->second.end(); i++) {
            add_clean_path_all_S(m->second, i, map_path_new);
        }
    }

    map_path=map_path_new;
    map_path.setName2NumQ(map_path_new.getName2NumQ());
    map_path.setName2NumS(map_path_new.getName2NumS());
    map_path.setNum2NameQ(map_path_new.getNum2NameQ());
    map_path.setNum2NameS(map_path_new.getNum2NameS());

    if (verbose > 0) {
        std::cout << "nb of matches: " << map_path.getNbMatchesInMapPath() << std::endl;
        std::cout << "nb of paths: " << map_path.getNbDistinctPaths() << std::endl;
    }
}
//---------------------------------------------------------------------------
void BLRMatchJoin::merge(BLRMatchPath &map_path, int verbose) {

    if (verbose > 0) {
        std::cout << "Merge on query (Compute score with length)." << std::endl;
    }

    BLRMatchPath map_path_merged;
    map_path_merged.setName2NumQ(map_path.getName2NumQ());
    map_path_merged.setName2NumS(map_path.getName2NumS());
    map_path_merged.setNum2NameQ(map_path.getNum2NameQ());
    map_path_merged.setNum2NameS(map_path.getNum2NameS());

    unsigned count_merge = 0;
    for (BLRMatchPath::MapPath::iterator m = map_path.begin(); m != map_path.end(); m++)
    {

        //Build overlap graph
        Graph<unsigned long> graph;


        std::map<unsigned long, std::list<RangePairSet>::const_iterator> idToRps;
        for (std::list<RangePairSet>::iterator lrp_it1 =  m->second.begin(); lrp_it1 != m->second.end(); lrp_it1++) {
            graph.add_node(lrp_it1->getId());
            lrp_it1->computeScoreWithLength();
            idToRps[lrp_it1->getId()] = lrp_it1;
            std::list<RangePairSet>::const_iterator lrp_it2 = lrp_it1;
            lrp_it2++;
            while (lrp_it2 != m->second.end()) {
                if (lrp_it1->overlapQ_length(*lrp_it2) >= merge_overlap) {
                    graph.add_edge(lrp_it1->getId(), lrp_it2->getId());
                }
                lrp_it2++;
            }
        }

        //Find connexe componante
        std::vector<std::vector<unsigned long> > vec;
        graph.connexComp(vec);

        // merge rps of in each connex comp
        std::list<RangePairSet> rpsListAfterMerge;
        for (std::vector<std::vector<unsigned long> >::iterator it_vec = vec.begin(); it_vec != vec.end(); it_vec++) {
            unsigned size = it_vec->size();

            // Get first rps from connex comp
            RangePairSet firstRps = *idToRps[(*it_vec)[0]];

            //Merge rps with other rps from same connex comp
            for (unsigned i = 1; i < size; i++) {
                RangePairSet rps = *(idToRps[(*it_vec)[i]]);
                firstRps.mergeQ(rps);
                firstRps.orientSubjects();
            }
            rpsListAfterMerge.push_back(firstRps);
        }


        //Insert new merged rps
        for (std::list<RangePairSet>::iterator it = rpsListAfterMerge.begin(); it != rpsListAfterMerge.end(); it++) {
            map_path_merged.insert(*it);
            count_merge++;
        }
    }
    map_path=map_path_merged;
    map_path.setName2NumQ(map_path_merged.getName2NumQ());
    map_path.setName2NumS(map_path_merged.getName2NumS());
    map_path.setNum2NameQ(map_path_merged.getNum2NameQ());
    map_path.setNum2NameS(map_path_merged.getNum2NameS());

    if (verbose > 0) {
        std::cout << "Merged results:" << count_merge << std::endl;
        std::cout << "nb of matches: " << map_path.getNbMatchesInMapPath() << std::endl;
        std::cout << "nb of paths: " << map_path.getNbDistinctPaths() << std::endl;
    }

}
//----------------------------------------------------------------------------
void BLRMatchJoin::add_split_path(std::list<RangePairSet> &rp_list, std::list<RangePairSet>::iterator iter, BLRMatchPath &map_path)
//add a rangePairSet and post process it removing split nest
{
    bool found_over = false;
    std::list<RangePairSet> lrp;
    lrp.push_back(*iter);

    for (BLRMatchPath::MapPath::iterator m = map_path.begin(); m != map_path.end(); m++) {
        for (std::list<RangePairSet>::iterator lrp_it = lrp.begin();
             lrp_it != lrp.end();
             lrp_it++)
            for (std::list<RangePairSet>::iterator iter_list = m->second.begin();
                 iter_list != m->second.end(); iter_list++)
                if (lrp_it->overlapQ(*iter_list)) {
                    // if overlap
                    if (lrp_it->getLength() >= 100 &&
                        lrp_it->inserted(*iter_list) &&
                        fabs(iter_list->getIdentity()
                             - lrp_it->getIdentity()) <= para.getIdTolerance())
                        // if inserted with no overlap, longer than 100bp, same identity --> do nothing
                        continue;
                    std::list<RangePairSet> lrp2;
                    if (lrp_it->split(*iter_list, lrp2)) {
                        // split overlaping rangePairSet
                        found_over = true;
                        for (std::list<RangePairSet>::iterator lrp2_it
                                = lrp2.begin(); lrp2_it != lrp2.end(); lrp2_it++)
                            if (lrp2_it->getRangeQ().getLength() > para.getLenFilter())
                                // Add splitted fragments
                                lrp.push_back(*lrp2_it);
                    }
                }
    }

    if (found_over) {
        for (std::list<RangePairSet>::iterator lrp_it = lrp.begin();
             lrp_it != lrp.end(); lrp_it++)
            if (!lrp_it->empty()
                && lrp_it->getRangeQ().getLength() > para.getLenFilter()) {
                // rangepairset has been modified by split() --> insert  in the right place (sorted by identity)
                std::list<RangePairSet>::iterator it
                        = std::lower_bound(iter, rp_list.end(),
                                           *lrp_it,
                                           RangePair::lessIdentity);
                if (it == iter)
                    it++;
                rp_list.insert(it, *lrp_it);
            }
    } else if (iter->getRangeQ().getLength() > para.getLenFilter())
        // No change on rangepairset, insert it normally by increasing coordinates
        map_path.insert(*iter);

}
//----------------------------------------------------------------------------
void BLRMatchJoin::split(BLRMatchPath &map_path, int verbose)
// split nested path
{

    BLRMatchPath map_path_new;
    map_path_new.setName2NumQ(map_path.getName2NumQ());
    map_path_new.setName2NumS(map_path.getName2NumS());
    map_path_new.setNum2NameQ(map_path.getNum2NameQ());
    map_path_new.setNum2NameS(map_path.getNum2NameS());


    for (BLRMatchPath::MapPath::iterator m = map_path.begin(); m != map_path.end(); m++) {
        m->second.sort(RangePair::lessIdentity);
        for (std::list<RangePairSet>::iterator i = m->second.begin();
             i != m->second.end(); i++) {
            add_split_path(m->second, i, map_path_new);
        }
    }

    map_path=map_path_new;
    map_path.setName2NumQ(map_path_new.getName2NumQ());
    map_path.setName2NumS(map_path_new.getName2NumS());
    map_path.setNum2NameQ(map_path_new.getNum2NameQ());
    map_path.setNum2NameS(map_path_new.getNum2NameS());


    if (verbose > 0) {
        std::cout << "nb of matches: " << map_path.getNbMatchesInMapPath() << std::endl;
        std::cout << "nb of paths: " << map_path.getNbDistinctPaths() << std::endl;
    }

}

