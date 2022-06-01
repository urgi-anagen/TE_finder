
/**
 *
 * BLRMatchMap.cpp
 *
 **/
#include <cstdlib>
#include <regex.h>
#include <fstream>
#include "FragAlign.h"
#include "BLRMatchMap.h"
#include <SDGError.h>
#include "BLRMatchMapLoader.h"

//----------------------------------------------------------------------------
void BLRMatchMap::insert(RangePair &rangePair) {
    //insert rangePair in the right place
    std::list<RangePair> &al_list
            = map_align[Key(rangePair.getRangeQ().getNumChr(),
                            rangePair.getRangeS().getNumChr())];

    std::list<RangePair>::iterator r = std::lower_bound(al_list.begin(),
                                                        al_list.end(), rangePair);
    if (rangePair.getRangeQ().getMin() != r->getRangeQ().getMin()
        || rangePair.getRangeQ().getMax() != r->getRangeQ().getMax()
        || rangePair.getRangeS().getMin() != r->getRangeS().getMin()
        || rangePair.getRangeS().getMax() != r->getRangeS().getMax()
        || rangePair.getE_value() != r->getE_value()
        || rangePair.getScore() != r->getScore()
        || rangePair.getIdentity() != r->getIdentity())
        al_list.insert(r, rangePair);
}
//----------------------------------------------------------------------------
void BLRMatchMap::readAlign(std::istringstream &streamName, int verbose)
//  read from a istringstream
{
    BLRMatchMapLoader blrmm = BLRMatchMapLoader();
    blrmm.readAlign(*this, streamName, verbose);
}
//----------------------------------------------------------------------------
//
void BLRMatchMap::loadAlign(SDGString filename, int verbose)
//  load from txt file
{
    BLRMatchMapLoader blrmm = BLRMatchMapLoader();
    blrmm.loadAlign(*this, filename, verbose);
}
//----------------------------------------------------------------------------
void BLRMatchMap::readPath(std::istringstream &streamName, int verbose) {
    BLRMatchMapLoader blrmm = BLRMatchMapLoader();
    blrmm.readPath(*this, streamName, verbose);
}
//---------------------------------------------------------------------------
//
void BLRMatchMap::loadPath(SDGString filename, int verbose) {
    BLRMatchMapLoader blrmm = BLRMatchMapLoader();
    blrmm.loadPath(*this, filename, verbose);
}
//----------------------------------------------------------------------------
void BLRMatchMap::add_clean(std::list<RangePair> &rp_list,
                            std::list<RangePair>::iterator iter)
//add a rangePair and post process it removing conflicting subjects
{
    // run_test_search_wSW for a conflicting subject
    bool found_over = false;
    std::list<RangePair> lrp; //list of modified and cleaned RangePair
    lrp.push_back(*iter);
    for (MapAlign::iterator m = map_align.begin(); m != map_align.end(); m++)
        if (m->first.first == iter->getRangeQ().getNumChr() &&
            m->first.second != iter->getRangeS().getNumChr()) {
            // check overlap only with a different subject

            for (std::list<RangePair>::iterator lrp_it = lrp.begin();
                 lrp_it != lrp.end();
                 lrp_it++) {
                for (std::list<RangePair>::iterator iter_list = m->second.begin();
                     iter_list != m->second.end();
                     iter_list++) {
                    if (lrp_it->getScore() < iter_list->getScore()
                        && lrp_it->overlapQ(*iter_list)) {
                        found_over = true;
                        RangePair rp = lrp_it->diffQ(*iter_list);
                        if (!rp.empty()
                            && rp.getRangeQ().getLength() > para.getLenFilter()) {
                            lrp.push_back(rp);
                        }
                    } //end if (...)
                } //end loop for
            }//end loop for
        } //end if

    if (found_over) // RangePair found to overlap (conflicts!)
    {
        for (std::list<RangePair>::iterator lrp_it = lrp.begin();
             lrp_it != lrp.end();
             lrp_it++)

            if (!lrp_it->empty()
                && lrp_it->getRangeQ().getLength() > para.getLenFilter()) {
                std::list<RangePair>::iterator it
                        = std::lower_bound(iter, rp_list.end(),
                                           *lrp_it,
                                           RangePair::greaterScore); // run_test_search_wSW for the right place to insert
                while (it != rp_list.end() && it == iter)
                    it++;
                rp_list.insert(it, *lrp_it);
            }
    } else // already cleaned RangePair
    if (!iter->empty() && iter->getRangeQ().getLength() > para.getLenFilter() && iter->getScore() > 0)
        insert(*iter);
}
//----------------------------------------------------------------------------
void BLRMatchMap::clean_conflicts(void)
// removing conflicting subjects
{
    std::list<RangePair> rp_list;
    for (MapAlign::iterator m = map_align.begin(); m != map_align.end(); m++) {
        while (!m->second.empty()) {
            RangePair rp = m->second.back();
            m->second.pop_back();
            if (para.getEvalFilter() >= rp.getE_value()
                || para.getIdFilter() <= rp.getIdentity()
                || para.getLenFilter() <= rp.getLength()) {
                rp_list.push_back(rp);
            }
        }
    }
    map_align.clear();

    rp_list.sort(RangePair::greaterScore);
    for (std::list<RangePair>::iterator i = rp_list.begin();
         i != rp_list.end(); i++)
        add_clean(rp_list, i);
}

//----------------------------------------------------------------------------
void BLRMatchMap::insert_path(RangePairSet &rangePair)
//insert a RangePairSet in map_path at the right place
{
    std::list<RangePairSet> &al_list
            = map_path[Key(rangePair.getRangeQ().getNumChr(),
                           rangePair.getRangeS().getNumChr())];

    std::list<RangePairSet>::iterator r = std::lower_bound(al_list.begin(),
                                                           al_list.end(), rangePair);
    if (rangePair.getRangeQ().getMin() != r->getRangeQ().getMin()
        || rangePair.getRangeQ().getMax() != r->getRangeQ().getMax()
        || rangePair.getRangeS().getMin() != r->getRangeS().getMin()
        || rangePair.getRangeS().getMax() != r->getRangeS().getMax()
        || rangePair.getE_value() != r->getE_value()
        || rangePair.getScore() != r->getScore()
        || rangePair.getIdentity() != r->getIdentity())
        al_list.insert(r, rangePair);
}

//----------------------------------------------------------------------------
void BLRMatchMap::add_clean_path_same_S(std::list<RangePairSet> &rp_list,
                                        std::list<RangePairSet>::iterator iter)
//add a range pair set and post process it removing conflicting subjects
{
    // run_test_search_wSW for a conflicting subject
    bool found_over = false;
    std::list<RangePairSet> &list
            = map_path[Key(iter->getRangeQ().getNumChr(),
                           iter->getRangeS().getNumChr())];
    for (std::list<RangePairSet>::iterator iter_list = list.begin();
         iter_list != list.end(); iter_list++)
        if (RangePair::greaterScore(*iter_list, *iter)
            && iter->overlapQ(*iter_list))
            if (iter->diffQ(*iter_list))
                found_over = true;

    if (found_over) {
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
    } else if (iter->getRangeQ().getLength() > para.getLenFilter()
               && iter->getScore() > 0)
        insert_path(*iter);

}
//----------------------------------------------------------------------------
void BLRMatchMap::add_clean_path_all_S(std::list<RangePairSet> &rp_list,
                                       std::list<RangePairSet>::iterator iter,
                                       int verbose)
//add a rangePairSet and post process it removing conflicting subjects
{
    bool found_over = false;
    bool atLeastOneOverlap = false;
    unsigned nbseqS = getNbSseq();
    //unsigned s=1;
    // start subject at -1 to consider merged data
    long s = -1;
    // for current query iter on each subject in map_path.
    while (s <= nbseqS) {
        // list of subject
        std::list<RangePairSet> &list
                = map_path[Key(iter->getRangeQ().getNumChr(), s)];
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
        insert_path(*iter);
    }
}
//----------------------------------------------------------------------------
void BLRMatchMap::clean_path(bool same_S, int verbose)
// removing conflicting subjects
{
    std::list<RangePairSet> rp_list;
    for (MapPath::iterator m = map_path.begin(); m != map_path.end(); m++) {
        while (!m->second.empty()) {
            RangePairSet rp = m->second.back();
            m->second.pop_back();
            rp.computeScoreWithLengthAndId();
            rp_list.push_back(rp);
        }
    }
    map_path.clear();
    rp_list.sort(RangePair::greaterScore);

    if (same_S)
        for (std::list<RangePairSet>::iterator i = rp_list.begin();
             i != rp_list.end(); i++)
            add_clean_path_same_S(rp_list, i);
    else
        for (std::list<RangePairSet>::iterator i = rp_list.begin();
             i != rp_list.end(); i++) {
            if (verbose > 0)
                std::cout << "Add " << numQ2name(i->getRangeQ().getNumChr()) << "-"
                          << numS2name(i->getRangeS().getNumChr())
                          << "-" << i->getScore() << std::endl << std::flush;
            add_clean_path_all_S(rp_list, i, verbose - 1);
            if (verbose > 0) {
                std::cout << "nb of matches: " << getNbMatchesInMapPath() << std::endl;
                std::cout << "nb of paths: " << getNbDistinctPaths() << std::endl;
            }
        }
}
//----------------------------------------------------------------------------
void BLRMatchMap::add_split_path(std::list<RangePairSet> &rp_list, std::list<RangePairSet>::iterator iter)
//add a rangePairSet and post process it removing split nest
{
    bool found_over = false;
    std::list<RangePairSet> lrp;
    lrp.push_back(*iter);

    for (MapPath::iterator m = map_path.begin(); m != map_path.end(); m++) {
        // Search same query chromosome
        if (m->first.first == iter->getRangeQ().getNumChr()) {
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
        } //end if
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
        insert_path(*iter);

}
//----------------------------------------------------------------------------
void BLRMatchMap::split_path(void)
// split nested path
{
    //transform map_path in rp_list
    std::list<RangePairSet> rp_list;
    for (MapPath::iterator m = path_begin(); m != path_end(); m++) {
        while (!m->second.empty()) {
            RangePairSet rp = m->second.back();
            m->second.pop_back();
            rp_list.push_back(rp);
        }
    }
    map_path.clear();
    // score no longer valid after split of RangePairSet by method
    // split
    rp_list.sort(RangePair::lessIdentity);
    for (std::list<RangePairSet>::iterator i = rp_list.begin();
         i != rp_list.end(); i++)
        add_split_path(rp_list, i);
}
/*//----------------------------------------------------------------------------
void BLRMatchMap::insert_path_static(MapPath mapPath, RangePairSet &rangePair)
//insert a RangePairSet in map_path at the right place
{
    std::list<RangePairSet> &al_list
            = mapPath[Key(rangePair.getRangeQ().getNumChr(),
                          rangePair.getRangeS().getNumChr())];

    std::list<RangePairSet>::iterator r = std::lower_bound(al_list.begin(),
                                                           al_list.end(), rangePair);
    if (rangePair.getRangeQ().getMin() != r->getRangeQ().getMin()
        || rangePair.getRangeQ().getMax() != r->getRangeQ().getMax()
        || rangePair.getRangeS().getMin() != r->getRangeS().getMin()
        || rangePair.getRangeS().getMax() != r->getRangeS().getMax()
        || rangePair.getE_value() != r->getE_value()
        || rangePair.getScore() != r->getScore()
        || rangePair.getIdentity() != r->getIdentity())
        al_list.insert(r, rangePair);
}*/
//----------------------------------------------------------------------------
void BLRMatchMap::mapPath(bool joining, bool clean_before, bool clean_after, bool merged, int verbose) {
    map_path.clear();

    if (joining) {
        if (verbose > 0)
            std::cout << "Join parameters: dist_pen=" << para.getDist_pen()
                      << " gap_pen=" << para.getGap_pen()
                      << " overlap=" << para.getOverlap() << std::endl << std::flush;


        for (MapAlign::iterator m = map_align.begin(); m != map_align.end(); m++) {
            FragAlign fragAlign(para.getDist_pen(), 0, para.getGap_pen(),
                                para.getOverlap());
            map_path[m->first] = fragAlign.join(m->second);
            m->second.clear();
        }
        if (verbose > 0) {
            std::cout << "After join:\nnb of matches: " << getNbMatchesInMapPath() << std::endl;
            std::cout << "nb of paths: " << getNbDistinctPaths() << std::endl;

            std::cout << "Write joined matches in BED format." << std::endl;
            SDGString filename = para.getPrefixFileName() + ".joined.bed";
            std::list<RangePairSet> copy_list = copyRpsListFromMapPath();
            SDGString color = "253,63,146";
            writeBED(filename, copy_list, verbose - 1);
        }
        if (merged) {
            if (verbose > 0)
                std::cout << "Compute score with length." << std::endl;
            computeScoreWithLength();

            // prepare rpsList for merge
       /*     unsigned long id = 1;
            rpsList.clear();
            for (MapPath::iterator m = path_begin(); m != path_end(); m++) {
                while (!m->second.empty()) {
                    RangePairSet rp = m->second.front();
                    m->second.pop_front();
                    rp.setId(id);
                    rpsList.push_back(rp);
                    id++;
                }
            }*/
            
            rpsList=getRpsListFromMapPath();
            map_path.clear();

            // merge
            if (verbose > 0)
                std::cout << "Merge on query." << std::endl;
            merge(verbose);
            if (verbose > 0) {
                std::cout << "Write merged matches in BED format." << std::endl;
                SDGString filename = para.getPrefixFileName() + ".merged.bed";
                std::list<RangePairSet> copy_list = copyRpsListFromMapPath();
                SDGString color = "0,255,0";
                writeBED(filename, copy_list, verbose - 1);
            }
        }
        if ((clean_before | clean_after) & (verbose > 0))
            std::cout << "Clean the connections..." << std::endl << std::flush;
        if (clean_before)
            clean_path(true, verbose - 1); // clean when same subject
        if (clean_after)
            clean_path(false, verbose - 1); // clean with all subjects
        if ((clean_before | clean_after) & (verbose > 0)) {
            std::cout << "After clean:\nnb of matches: " << getNbMatchesInMapPath() << std::endl;
            std::cout << "nb of paths: " << getNbDistinctPaths() << std::endl;
            std::cout << "Connections were cleaned when necessary." << std::endl;

            std::cout << "Write cleaned matches in BED format." << std::endl;
            SDGString filename = para.getPrefixFileName() + ".cleaned.bed";
            std::list<RangePairSet> copy_list = copyRpsListFromMapPath();
            SDGString color = "255,0,0";
            writeBED(filename, copy_list, verbose - 1);
        }
        if (clean_before || clean_after) {
            if (verbose > 0)
                std::cout << "Split the connections..." << std::endl << std::flush;
            split_path();
            if (verbose > 0) {
                std::cout << "nb of matches: " << getNbMatchesInMapPath() << std::endl;
                std::cout << "nb of paths: " << getNbDistinctPaths() << std::endl;
                std::cout << "Connections were splitted when necessary." << std::endl;
                std::cout << "Write split matches in BED format." << std::endl;
                SDGString filename = para.getPrefixFileName() + ".split.bed";
                std::list<RangePairSet> copy_list = copyRpsListFromMapPath();
                SDGString color = "237,127,16";
                writeBED(filename, copy_list, verbose - 1);
            }
        }
    } else  { // no join
        int count = 0;
        if (verbose > 0)
            std::cout << "No join, considering ";
        for (MapAlign::iterator m = map_align.begin(); m != map_align.end(); m++) {
            std::list<RangePairSet> path;
            for (std::list<RangePair>::iterator i = m->second.begin();
                 i != m->second.end(); i++) {
                count++;
                path.push_back(RangePairSet(*i));
            }
            map_path[m->first] = path;
            m->second.clear();
        }
        if (verbose > 0)
            std::cout << count << " matches" << std::endl;
    }
    map_align.clear();
}
//----------------------------------------------------------------------------
void BLRMatchMap::mapPathJoinOnlyForTest(bool joining, bool clean_before, bool clean_after, int verbose) {
    map_path.clear();

    if (joining) {
        if (verbose > 0)
            std::cout << "Join parameters: dist_pen=" << para.getDist_pen()
                      << " gap_pen=" << para.getGap_pen()
                      << " overlap=" << para.getOverlap() << std::endl << std::flush;
        for (MapAlign::iterator m = map_align.begin(); m != map_align.end(); m++) {
            FragAlign fragAlign(para.getDist_pen(), 0, para.getGap_pen(),
                                para.getOverlap());
            map_path[m->first] = fragAlign.join(m->second);
            m->second.clear();
        }
        if (verbose > 0) {
            std::cout << "nb of matches: " << getNbMatchesInMapPath() << std::endl;
            std::cout << "nb of paths: " << getNbDistinctPaths() << std::endl;
        }
    }
    return;
}

//----------------------------------------------------------------------------
std::list<RangePairSet> BLRMatchMap::getRpsList(void) // Build list of RangePairSet
{
    std::list<RangePairSet> rps_list;
    for (std::list<RangePairSet>::iterator it = rpsList.begin(); it != rpsList.end(); it++) {
        RangePairSet rps = *it;
        std::string subject_name;
        if (rps.getNumSubject() == -1)
            subject_name = "-1";
        else
            subject_name = num2nameS[rps.getNumSubject()];
        rps.setQSName(num2nameQ[rps.getNumQuery()],subject_name, num2nameS );
        rps_list.push_back(rps);
    }
    return rps_list;
};

//----------------------------------------------------------------------------
// Note: This methode is destructive for mapPath
std::list<RangePairSet> BLRMatchMap::getRpsListFromMapPath(void) // Build list of RangePairSet
{
    unsigned count=0;
    std::list<RangePairSet> rps_list;
    for (MapPath::iterator m = path_begin(); m != path_end(); m++) {
        while (!m->second.empty()) {
            RangePairSet rp = m->second.front();
            rp.setId(++count);
            m->second.pop_front();
            rps_list.push_back(rp);
        }
    }
    return rps_list;
}
//----------------------------------------------------------------------------
std::list<RangePairSet> BLRMatchMap::copyRpsListFromMapPath(void) {
    std::list<RangePairSet> copy_rps_list;
    for (MapPath::iterator m = path_begin(); m != path_end(); m++) {
        for (std::list<RangePairSet>::iterator it = m->second.begin(); it != m->second.end(); it++) {
            copy_rps_list.push_back(*it);
        }
    }
    return copy_rps_list;
}
/*//---------------------------------------------------------------------------
void BLRMatchMap::selectQregex(SDGString regex) {
    regex_t preg;
    regcomp(&preg, regex, REG_EXTENDED | REG_NEWLINE | REG_NOSUB);

    MapAlign::iterator iter_hash, prev;  // iterator to visit hash_align
    iter_hash = begin();
    while (iter_hash != end()) {
        prev = iter_hash++;
        if (regexec(&preg, SDGString(num2nameQ[prev->first.first]), 1, NULL, 0) != 0) {
            map_align.erase(prev);
        }
    }
    regfree(&preg);
}

//---------------------------------------------------------------------------
void BLRMatchMap::select(bool subject, bool clean_before, bool clean_after) {
    std::ostringstream outfile;
    if (subject) {
        RangeMap matchmap;
        MapPath::iterator iter_hash;  // iterator to visit hash_align
        std::list<RangePairSet>::iterator iter_list;  //iterator to visit align_list

        iter_hash = path_begin();
        while (iter_hash != path_end()) {
            iter_list = iter_hash->second.begin();
            while (iter_list != iter_hash->second.end()) {
                std::string chrS_name = num2nameS[(iter_list->getRangeS().getNumChr()) - 1];
                matchmap.add(RangeSeq("", chrS_name, 0, 0));

                iter_list++;
            }
            iter_hash++;
        }

        if (clean_before || clean_after)
            outfile << para.getParameterFileName().beforelast(".param")
                    << ".clean_match.subject_selected";
        else
            outfile << para.getParameterFileName().beforelast(".param")
                    << ".match.subject_selected";
        matchmap.selectSrcSeq(outfile.str().c_str(), subject_db);
    } else {

        RangeMap matchmap;
        MapPath::iterator iter_hash;  // iterator to visit hash_align
        std::list<RangePairSet>::iterator iter_list;  //iterator to visit align_list
        iter_hash = path_begin();
        while (iter_hash != path_end()) {
            iter_list = iter_hash->second.begin();
            while (iter_list != iter_hash->second.end()) {
                std::string chrQ_name = num2nameQ[(iter_list->getRangeQ().getNumChr()) - 1];
                matchmap.add(RangeSeq("", chrQ_name, 0, 0));

                iter_list++;
            }
            iter_hash++;
        }

        if (clean_before || clean_after)
            outfile << para.getParameterFileName().beforelast(".param")
                    << ".clean_match.query_selected";
        else
            outfile << para.getParameterFileName().beforelast(".param")
                    << ".match.query_selected";
        matchmap.selectSrcSeq(outfile.str().c_str(), query_db);
    }
    std::cout << "ok!" << std::endl;
}*/

//---------------------------------------------------------------------------
void BLRMatchMap::writeMatch(const SDGString &filename, int verbose) {
    if (verbose > 0)
        std::cout << "writing 'tab' file..." << std::flush;
    std::ofstream fout(filename);

    fout << "query.name"
         << "\t" << "query.start"
         << "\t" << "query.end"
         << "\t" << "query.length"
         << "\t" << "query.length.%"
         << "\t" << "match.length.%"
         << "\t" << "subject.name"
         << "\t" << "subject.start"
         << "\t" << "subject.end"
         << "\t" << "subject.length"
         << "\t" << "subject.length.%"
         << "\t" << "E.value"
         << "\t" << "Score"
         << "\t" << "Identity"
         << "\t" << "path"
         << std::endl;

    unsigned path_id = 0;
    for (MapPath::iterator iter_hash = path_begin();
         iter_hash != path_end(); iter_hash++) {

        std::string query_name = num2nameQ[iter_hash->first.first];
        std::string subject_name;
        if (same_db)
            subject_name = num2nameQ[iter_hash->first.second];
        else
            subject_name = num2nameS[iter_hash->first.second];

        std::string subseqname(subject_name, 0, subject_name.find(" "));
        unsigned querylen = query_db.find(query_name).length();
        unsigned subjectlen = subject_db.find(subject_name).length();

        for (std::list<RangePairSet>::iterator iter_list
                = iter_hash->second.begin(); iter_list != iter_hash->second.end();
             iter_list++) {
            //iter_list->view();

            RangeAlignSet rasQ = iter_list->getRangeAlignSetQ();
            RangeAlignSet rasS = iter_list->getRangeAlignSetS();


            unsigned rangeQlen = rasQ.getLengthSet();
            unsigned rangeSlen = rasS.getLengthSet();

            fout << query_name
                 << "\t" << rasQ.getStart()
                 << "\t" << rasQ.getEnd()
                 << "\t" << rangeQlen
                 << "\t" << (double) (rangeQlen) / querylen
                 << "\t" << (double) (rangeQlen) / subjectlen
                 << "\t" << subseqname
                 << "\t" << rasS.getStart()
                 << "\t" << rasS.getEnd()
                 << "\t" << rangeSlen
                 << "\t" << (double) (rangeSlen) / subjectlen
                 << "\t" << iter_list->getE_value()
                 << "\t" << iter_list->getScore()
                 << "\t" << iter_list->getIdentity()
                 << "\t" << ++path_id
                 << std::endl;
        }
    }
    if (verbose > 0)
        std::cout << " done" << std::endl;
}
//---------------------------------------------------------------------------
void BLRMatchMap::writePath(std::ostream& out, int verbose) {
    if (verbose > 0)
        std::cout << "writing 'path' file..." << std::flush;

    unsigned path_id = 0;
    for (MapPath::iterator m = path_begin(); m != path_end(); m++) {
        for (std::list<RangePairSet>::iterator it = m->second.begin(); it != m->second.end(); it++) {
            std::string query_name = num2nameQ[it->getNumQuery()];
            unsigned id = ++path_id;
            it->write(out, id, query_name, num2nameS);
        }
    }

    if (verbose > 0)
        std::cout << " done" << std::endl;
}
//---------------------------------------------------------------------------
void BLRMatchMap::writePathAttr(std::ostream& out, int verbose) {
    if (verbose > 0)
        std::cout << "writing 'path' file..." << std::flush;

    unsigned path_id = 0;
    for (MapPath::iterator m = path_begin(); m != path_end(); m++) {
        for (std::list<RangePairSet>::iterator it = m->second.begin(); it != m->second.end(); it++) {
            std::string query_name = num2nameQ[it->getNumQuery()];
            unsigned id = ++path_id;
            it->writeRpsAttr(out, id, query_name, num2nameS);
        }
    }

    if (verbose > 0)
        std::cout << " done" << std::endl;
}
//---------------------------------------------------------------------------
void BLRMatchMap::writePath(const SDGString &filename,
                            std::list<RangePairSet> &rps_list, int verbose) {
    if (verbose > 0)
        std::cout << "writing 'path' file..." << std::flush;

    std::ofstream fout(filename);
    std::ofstream foutRpsAttr(filename + ".attr");

    writeRpsList(rps_list, fout);
    writeRpsListAttribute(rps_list, foutRpsAttr);

    if (verbose > 0)
        std::cout << " done" << std::endl;
}

//---------------------------------------------------------------------------
void BLRMatchMap::writeRpsListAttribute(std::list<RangePairSet> &rps_list, std::ostream &out) {
    unsigned path_id = 0;
    for (std::list<RangePairSet>::iterator iter_list
            = rps_list.begin(); iter_list != rps_list.end();
         iter_list++) {
        std::string query_name = num2nameQ[iter_list->getNumQuery()];
        unsigned id = ++path_id;
        iter_list->writeRpsAttr(out, id, query_name, num2nameS);
    }
}

//---------------------------------------------------------------------------
void BLRMatchMap::writeRpsList(std::list<RangePairSet> &rps_list, std::ostream &out) {
    unsigned path_id = 0;
    for (std::list<RangePairSet>::iterator iter_list
            = rps_list.begin(); iter_list != rps_list.end();
         iter_list++) {
        std::string query_name = num2nameQ[iter_list->getNumQuery()];
        unsigned id = ++path_id;
        iter_list->write(out, id, query_name, num2nameS);
    }
}

//---------------------------------------------------------------------------
void BLRMatchMap::writeBED(const SDGString &filename, const std::list<RangePairSet> &rps_list,
                           int verbose) {
    std::ostringstream bedStream;
    writeBED(bedStream, rps_list, verbose);
    std::ofstream bedFile(filename);
    bedFile << bedStream.str();
}

//---------------------------------------------------------------------------
void
BLRMatchMap::writeBED(std::ostream &out, const std::list<RangePairSet> &rps_list, int verbose) {
    if (verbose > 0)
        std::cout << "writing 'bed' file..." << std::flush;

    if (verbose > 0) {
        std::cout << " " << std::endl;
        std::cout << "writeBED rpsList size " << rps_list.size() << std::endl;
    }


    for (std::list<RangePairSet>::const_iterator iter_list
            = rps_list.begin(); iter_list != rps_list.end();
         iter_list++) {
        std::string query_name = num2nameQ[iter_list->getNumQuery()];
        iter_list->writeBED(out, query_name, num2nameS);
    }
}
//---------------------------------------------------------------------------
void BLRMatchMap::writeGFF3(const SDGString &filename, const std::list<RangePairSet> &rps_list,
                           int verbose) {
    std::ostringstream gffStream;
    writeGFF3(gffStream, rps_list, verbose);
    std::ofstream bedFile(filename);
    bedFile << gffStream.str();
}

//---------------------------------------------------------------------------
void
BLRMatchMap::writeGFF3(std::ostream &out, const std::list<RangePairSet> &rps_list, int verbose) {
    if (verbose > 0)
        std::cout << "writing 'gff3' file..." << std::flush;

    if (verbose > 0) {
        std::cout << " " << std::endl;
        std::cout << "writeGFF3 rpsList size " << rps_list.size() << std::endl;
    }

    unsigned id=0;
    for (std::list<RangePairSet>::const_iterator iter_list
            = rps_list.begin(); iter_list != rps_list.end();
         iter_list++) {
        std::string query_name = num2nameQ[iter_list->getNumQuery()];
        iter_list->writeGFF3(out, ++id, query_name, num2nameS);
    }
}
/*

//---------------------------------------------------------------------------
void BLRMatchMap::writePathForMergedS(std::ostream &fout, unsigned path_id, std::string query_name,
                                      std::list<RangePair> path) {
    for (std::list<RangePair>::iterator i = path.begin(); i != path.end(); i++) {

        long numChrS = i->getRangeS().getNumChr();

        std::string nameS = num2nameS[numChrS];
        fout << path_id << "\t" << query_name
             << "\t" << i->getRangeQ().getStart()
             << "\t" << i->getRangeQ().getEnd()
             << "\t" << nameS
             << "\t" << i->getRangeS().getStart()
             << "\t" << i->getRangeS().getEnd()
             << "\t" << i->getE_value()
             << "\t" << i->getScore()
             << "\t" << i->getIdentity()
             << std::endl;
    }
}
*/

//---------------------------------------------------------------------------
RangeMap BLRMatchMap::writeMap(const SDGString &filename, int verbose) {
    if (verbose > 0)
        std::cout << "writing 'map' file..." << std::flush;

    RangeMap matchmap;
    unsigned path_id = 0;
    for (MapPath::iterator iter_hash = path_begin();
         iter_hash != path_end(); iter_hash++) {
        std::string query_name = num2nameQ[iter_hash->first.first];
        std::string subject_name;
        if (same_db)
            subject_name = num2nameQ[iter_hash->first.second];
        else
            subject_name = num2nameS[iter_hash->first.second];

        std::string subseqname(subject_name, 0, subject_name.find(" "));


        for (std::list<RangePairSet>::iterator iter_list
                = iter_hash->second.begin(); iter_list != iter_hash->second.end();
             iter_list++) {
            RangeAlignSet rasQ = iter_list->getRangeAlignSetQ();
            RangeAlignSet rasS = iter_list->getRangeAlignSetS();

            SDGString copyname = subseqname + "." + SDGString(++path_id);

            if (!rasS.isPlusStrand())
                rasQ.reverse();

            matchmap.add(RangeSeq(rasQ, copyname, query_name));
        }
    }

    matchmap.save(filename);
    if (verbose > 0)
        std::cout << " done" << std::endl;
    return matchmap;
}

//---------------------------------------------------------------------------
void BLRMatchMap::writeSeq(const RangeMap &matchmap, const SDGString &filename, int verbose) {
    if (verbose > 0)
        std::cout << "writing 'fasta' file..." << std::endl;
    matchmap.writeSeq(filename, para.getQuery());
    if (verbose > 0)
        std::cout << " done" << std::endl;
}

//---------------------------------------------------------------------------
void BLRMatchMap::writeMapAlign(std::ostream &out) {
    for (MapAlign::iterator m = map_align.begin(); m != map_align.end(); m++) {
        while (!m->second.empty()) {
            Key key = m->first;
            SDGString numQuery = SDGString(key.first);
            SDGString numSubject = SDGString(key.second);
            //std::cout<<" "<<std::endl;
            //std::cout<<"Num query: "<<numQuery<<" Num subject: "<<numSubject<<std::endl;
            RangePair range_pair = m->second.back();
            SDGString numChr = SDGString(range_pair.getRangeQ().getNumChr());
            SDGString nameSeq = range_pair.getRangeQ().getNameSeq();
            //std::cout<<"Query Name Seq: "<<nameSeq<<"Query Num chr: "<<numChr<<std::endl;
            m->second.pop_back();
            range_pair.write(out);
        }
    }

}

//---------------------------------------------------------------------------
void BLRMatchMap::contigOverlap(void) {
    if (para.getQuery() != para.getBank()) {
        std::cout << "MATCHER " << para.getQuery() << " vs " << para.getBank()
                  << " -> no contig overlap !!" << std::endl;
        return;
    }


    unsigned count = 0;
    RangeMap matchmap;
    MapPath::iterator iter_hash = map_path.begin();
    while (iter_hash != map_path.end()) {
        std::list<RangePairSet>::iterator iter_list = iter_hash->second.begin();
        while (iter_list != iter_hash->second.end()) {
            if (
                    (
                            iter_list->getRangeQ().getMin() == 1
                            && iter_list->getRangeS().getMax()
                               == query_db[(iter_list->getRangeS().getNumChr()) - 1].length()
                    )
                    ||
                    (
                            iter_list->getRangeS().getMin() == 1
                            && iter_list->getRangeQ().getMax()
                               == query_db[(iter_list->getRangeQ().getNumChr()) - 1].length()
                    )
                    ) {
                std::cout << *iter_list << std::endl;

                SDGString name(++count);
                std::string chr_name = num2nameQ[(iter_list->getRangeQ().getNumChr()) - 1];

                if (iter_list->getRangeQ().isPlusStrand())
                    matchmap.add(RangeSeq(name, chr_name,
                                          iter_list->getRangeQ().getStart(),
                                          iter_list->getRangeQ().getEnd()));
                else
                    matchmap.add(RangeSeq(name, chr_name,
                                          iter_list->getRangeQ().getEnd(),
                                          iter_list->getRangeQ().getStart()));

                chr_name = num2nameS[(iter_list->getRangeS().getNumChr()) - 1];

                if (iter_list->getRangeS().isPlusStrand())
                    matchmap.add(RangeSeq(name, chr_name,
                                          iter_list->getRangeS().getStart(),
                                          iter_list->getRangeS().getEnd()));
                else
                    matchmap.add(RangeSeq(name, chr_name,
                                          iter_list->getRangeS().getEnd(),
                                          iter_list->getRangeS().getStart()));
            }

            iter_list++;
        }
        iter_hash++;
    }
}

//---------------------------------------------------------------------------
unsigned BLRMatchMap::getNbSseq(void) {
    if (same_db)
        return num2nameQ.size();
    else
        return num2nameS.size();
}
//---------------------------------------------------------------------------
unsigned BLRMatchMap::getNbMatchesInMapAlign(void) {
    unsigned nbMatches = 0;
    for (MapAlign::iterator m = map_align.begin(); m != map_align.end(); m++)
        nbMatches += m->second.size();
    return nbMatches;
}

//---------------------------------------------------------------------------
unsigned BLRMatchMap::getNbMatchesInMapPath(void) {
    unsigned nbMatches = 0;
    for (MapPath::iterator m = map_path.begin(); m != map_path.end(); m++)
        for (std::list<RangePairSet>::iterator i = m->second.begin(); i != m->second.end(); i++)
            nbMatches += i->getNbRangePairs();
    return nbMatches;
}

//---------------------------------------------------------------------------
unsigned BLRMatchMap::getNbDistinctPaths(void) {
    unsigned nbPaths = 0;
    for (MapPath::iterator m = map_path.begin(); m != map_path.end(); m++)
        nbPaths += m->second.size();
    return nbPaths;
}
//---------------------------------------------------------------------------
void BLRMatchMap::computeScoreWithLength(std::list<RangePairSet> &rpsList) {
    for (std::list<RangePairSet>::iterator i = rpsList.begin(); i != rpsList.end(); i++) {
        i->computeScoreWithLengthAndId();
    }
}
//---------------------------------------------------------------------------
void BLRMatchMap::computeScoreWithLength() {
    for (MapPath::iterator m = map_path.begin(); m != map_path.end(); m++) {
        computeScoreWithLength(m->second);
    }
}
//---------------------------------------------------------------------------
void BLRMatchMap::merge(int verbose) {

    // Assume that mapPath content is set in rpsList
    Graph<unsigned long> graph;


    std::map<unsigned long, std::list<RangePairSet>::const_iterator> idToRps;
    for (std::list<RangePairSet>::const_iterator lrp_it1 = rpsList.begin(); lrp_it1 != rpsList.end(); lrp_it1++) {
        RangePairSet rps1 = *lrp_it1;
        idToRps[lrp_it1->getId()]=lrp_it1;
        graph.add_node(lrp_it1->getId());
        std::list<RangePairSet>::const_iterator lrp_it2 = lrp_it1;
        lrp_it2++;
        while(lrp_it2 != rpsList.end()){
            RangePairSet rps2 = *lrp_it2;
            if (rps1.overlapQ_length(rps2) >= merge_overlap) {
                graph.add_edge(rps1.getId(), rps2.getId());
            }
            lrp_it2++;
        }
    }

    if(verbose>2) graph.view();
    //Merge according to connexe componant
    if(verbose>0){std::cout<<"run_test_search_wSW connexe componants"<<std::endl;}
    std::vector<std::vector<unsigned long> > vec;
    graph.connexComp(vec);
    if(verbose>0){std::cout<<"Nb connexe componante:"<<vec.size()<<std::endl;}
    if(verbose>0){std::cout<<"merge connexe componants"<<std::endl;}

    // merge rps of in each connex comp
    std::list<RangePairSet> rpsListAfterMerge;
    for (std::vector<std::vector<unsigned long> >::iterator it_vec = vec.begin(); it_vec != vec.end(); it_vec++) {
        unsigned size = it_vec->size();
        if(verbose>1) std::cout<<"** new connexe componant\n+ node "<<(*it_vec)[0]<<std::endl;
        // Get first rps from connex comp
        RangePairSet firstRps = *idToRps[(*it_vec)[0]];
        if(verbose>1) {std::cout<<"-first path : "<<std::flush; firstRps.view();};
        //Merge rps with other rps from same connex comp
        for (unsigned i = 1; i < size; i++) {
            if(verbose>1) std::cout<<"+ node "<<(*it_vec)[i]<<std::endl;
            RangePairSet rps=*(idToRps[(*it_vec)[i]]);
            if(verbose>1) {std::cout<<"-new path : "<<std::flush; rps.view();};
            firstRps.mergeQ(rps);
            firstRps.orientSubjects();
            if(verbose>1) {std::cout<<"-merged path : "<<std::flush; firstRps.view();};
        }
        rpsListAfterMerge.push_back(firstRps);
    }

    int count=0;
    //Insert new merged rps
    if(verbose>0){std::cout<<"store paths"<<std::endl;}
    for (std::list<RangePairSet>::iterator it = rpsListAfterMerge.begin(); it != rpsListAfterMerge.end(); it++) {
        insert_path(*it);
        count++;
    }

    if(verbose>0){std::cout<<"Nb merged paths:"<<count<<std::endl;}
}

//---------------------------------------------------------------------------
void BLRMatchMap::insert(RangePairSet &rp) {
    unsigned long currentId, previousId = 0;
    RangePairSet rps;
    std::list<RangePair> path;

    if (rpsList.size() != 0) {
        rps = rpsList.back();
        previousId = rps.getId();
    }

    currentId = rp.getId();
    if (currentId == previousId) {
        rps = rpsList.back();
        rpsList.pop_back();
        path = rps.getPath();
    }

    path.push_back(rp);
    rps.setRpsFromRpList(path);
    rpsList.push_back(rps);

}
/*

//----------------------------------------------------------------------------
void
BLRMatchMap::mapPathJoinAndComputeScoreWithLengthOnly(bool joining, bool clean_before, bool clean_after, int verbose) {
    map_path.clear();

    if (joining) {
        if (verbose > 0)
            std::cout << "Join parameters: dist_pen=" << para.getDist_pen()
                      << " gap_pen=" << para.getGap_pen()
                      << " overlap=" << para.getOverlap() << std::endl << std::flush;
        for (MapAlign::iterator m = map_align.begin(); m != map_align.end(); m++) {
            FragAlign fragAlign(para.getDist_pen(), 0, para.getGap_pen(),
                                para.getOverlap());
            map_path[m->first] = fragAlign.join(m->second);
            m->second.clear();
        }
        if (verbose > 0) {
            std::cout << "nb of matches: " << getNbMatchesInMapPath() << std::endl;
            std::cout << "nb of paths: " << getNbDistinctPaths() << std::endl;
        }

        if ((clean_before | clean_after) & (verbose > 0))
            std::cout << "Clean the connections..." << std::endl << std::flush;
        if (clean_before)
            clean_path(true, verbose - 1); // clean when same subject
        if (clean_after)
            clean_path(false, verbose - 1); // clean with all subjects
        if ((clean_before | clean_after) & (verbose > 0)) {
            std::cout << "nb of matches: " << getNbMatchesInMapPath() << std::endl;
            std::cout << "nb of paths: " << getNbDistinctPaths() << std::endl;
            std::cout << "Connections were cleaned when necessary." << std::endl;
        }
        if (verbose > 0) {
            std::cout << "Recompute score with length... " << getNbMatchesInMapPath() << std::endl;
        }

        for (MapPath::iterator m = map_path.begin(); m != map_path.end(); m++) {
            computeScoreWithLengthAndId(m->second);
        }
    }
    return;
}
*/
