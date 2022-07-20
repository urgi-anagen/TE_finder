//
// Created by Hadi Quesneville on 2018-12-21.
//

#include "BLRMatchAlign.h"

//---------------------------------------------------------------------------
void BLRMatchAlign::read(const BLRJoinParameter& param, std::istream& input_align, int verbose)
{
    unsigned countseqS=0,countseqQ=0;
    para=param;

    getName2NumQ().clear();
    getName2NumS().clear();
    getNum2NameQ().clear();
    getNum2NameS().clear();


    //Check format
    std::string str;
    char buff[1024];
    input_align.getline(buff,1023,'\n');
    str=buff;
    size_t n = std::count(str.begin(), str.end(), '\t');
    if (n!=8)
    {
        std::cout<<"BLRMatchAlign ERROR: Number of columns is:"<<n+1<<std::endl
                 <<"first line:\n"<<str<<std::endl
                 <<"! not an *align* formated file"<<std::endl;
        exit(0);
    }

    //Read the first used to test the file format
    std::istringstream is(str);
    RangePair rp;
    rp.readtxt(is);
    if(para.getEvalFilter()>rp.getE_value()
       && para.getIdFilter()<rp.getIdentity()
       && para.getLenFilter()<rp.getLength())
    {
        std::map<std::string,long>::iterator it
                =name2numQ.find(rp.getRangeQ().getNameSeq());
        if(it==name2numQ.end())
        {
            name2numQ[rp.getRangeQ().getNameSeq()]=++countseqQ;
            num2nameQ[countseqQ]=rp.getRangeQ().getNameSeq();
            rp.getRangeQ().setNumChr(countseqQ);
        }

        rp.getRangeQ().setNumChr(it->second);

        it=name2numS.find(rp.getRangeS().getNameSeq());
        if(it==name2numS.end())
        {
            name2numS[rp.getRangeS().getNameSeq()]=++countseqS;
            num2nameS[countseqS]=rp.getRangeS().getNameSeq();
            rp.getRangeS().setNumChr(countseqS);
        }
        else
            rp.getRangeS().setNumChr(it->second);
        insert(rp);
    }

    // read the others
    while(input_align)
    {
        RangePair rp;
        rp.readtxt(input_align);
        if(input_align)
        {
            if(para.getEvalFilter()<rp.getE_value()
               || para.getIdFilter()>rp.getIdentity()
               || para.getLenFilter()>rp.getLength())
                continue;

            std::map<std::string,long>::iterator it
                    =name2numQ.find(rp.getRangeQ().getNameSeq());
            if(it==name2numQ.end())
            {
                name2numQ[rp.getRangeQ().getNameSeq()]=++countseqQ;
                num2nameQ[countseqQ]=rp.getRangeQ().getNameSeq();
                rp.getRangeQ().setNumChr(countseqQ);
            }
            else
                rp.getRangeQ().setNumChr(it->second);


            it=name2numS.find(rp.getRangeS().getNameSeq());
            if(it==name2numS.end())
            {
                name2numS[rp.getRangeS().getNameSeq()]=++countseqS;
                num2nameS[countseqS]=rp.getRangeS().getNameSeq();
                rp.getRangeS().setNumChr(countseqS);
            }
            else
                rp.getRangeS().setNumChr(it->second);

            insert(rp);
        }
    }


    if(verbose>0)
    {
        std::cout<<"nb of matches: "<<getNbMatchesInMapAlign()<<std::endl;
        std::cout<<"nb of distinct queries: "<<getNbQseq()<<std::endl;
        std::cout<<"nb of distinct subjects: "<<getNbSseq()<<std::endl;
    }
}
//----------------------------------------------------------------------------
void BLRMatchAlign::setFromRpsList(const BLRJoinParameter& param, const std::list<RangePair>& rp_list, int verbose)
{
    unsigned countseqS=0,countseqQ=0;

    getName2NumQ().clear();
    getName2NumS().clear();
    getNum2NameQ().clear();
    getNum2NameS().clear();


    // read the others
    for(std::list<RangePair>::const_iterator rp_it=rp_list.begin(); rp_it != rp_list.end(); rp_it++)
    {

            RangePair rp(*rp_it);
            std::map<std::string,long>::iterator it
                    =name2numQ.find(rp.getRangeQ().getNameSeq());
            if(it==name2numQ.end())
            {
                name2numQ[rp.getRangeQ().getNameSeq()]=++countseqQ;
                num2nameQ[countseqQ]=rp.getRangeQ().getNameSeq();
                rp.getRangeQ().setNumChr(countseqQ);
            }
            else
            rp.getRangeQ().setNumChr(it->second);


            it=name2numS.find(rp.getRangeS().getNameSeq());
            if(it==name2numS.end())
            {
                name2numS[rp.getRangeS().getNameSeq()]=++countseqS;
                num2nameS[countseqS]=rp.getRangeS().getNameSeq();
                rp.getRangeS().setNumChr(countseqS);
            }
            else
                rp.getRangeS().setNumChr(it->second);

            insert(rp);
    }

    if(verbose>0)
    {
        std::cout<<"nb of matches: "<<getNbMatchesInMapAlign()<<std::endl;
        std::cout<<"nb of distinct queries: "<<getNbQseq()<<std::endl;
        std::cout<<"nb of distinct subjects: "<<getNbSseq()<<std::endl;
    }
}
//----------------------------------------------------------------------------
void BLRMatchAlign::insert(RangePair &rangePair) {
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
void BLRMatchAlign::add_clean_overlap(std::list<RangePair> &rp_list,
                                      std::list<RangePair>::iterator iter)
//add a rangePair and post process it removing conflicting subjects
{
    // search for a conflicting subject
    bool found_over = false;
    std::list<RangePair> lrp; //list of modified and cleaned RangePair
    lrp.push_back(*iter);
    for (MapAlign::iterator m = map_align.begin(); m != map_align.end(); m++) // strange loop !!!
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
                                           RangePair::greaterScore); // search for the right place to insert
                while (it != rp_list.end() && it == iter)
                    it++;
                rp_list.insert(it, *lrp_it);
            }
    } else // already cleaned RangePair
    if (!iter->empty() && iter->getRangeQ().getLength() > para.getLenFilter() && iter->getScore() > 0)
        insert(*iter);
}
//----------------------------------------------------------------------------
void BLRMatchAlign::add_clean_included(std::list<RangePair>::iterator iter)
//add a rangePair and post process it removing conflicting subjects
{
    // search for a conflicting subject
    bool found_over = false;
    for (MapAlign::iterator m = map_align.begin(); m != map_align.end(); m++)
        if (m->first.first == iter->getRangeQ().getNumChr() &&
            m->first.second != iter->getRangeS().getNumChr()) {
            for (std::list<RangePair>::iterator iter_list = m->second.begin();
                 iter_list != m->second.end();
                 iter_list++) {
                if (iter->getRangeS().getNumChr() != iter_list->getRangeS().getNumChr()
                    && iter->getScore() < iter_list->getScore()
                    && iter_list->includedQ(*iter)) {
                    found_over = true;
                    break;
                } //end if (...)
            } //end loop for
        } //end if


    if (!found_over) // RangePair not found to overlap (no conflicts!)
        insert(*iter);
}
//----------------------------------------------------------------------------
void BLRMatchAlign::clean_conflicts(void)
// removing conflicting subjects
{
    std::vector< std::list<RangePair> > vec_rp_list(getNbQseq());
    for (MapAlign::iterator m = map_align.begin(); m != map_align.end(); m++) {
        while (!m->second.empty()) {
            RangePair rp = m->second.back();
            m->second.pop_back();
            if (para.getEvalFilter() >= rp.getE_value()
                || para.getIdFilter() <= rp.getIdentity()
                || para.getLenFilter() <= rp.getLength()) {
                vec_rp_list[m->first.first-1].push_back(rp);
            }
        }
    }
    map_align.clear();

    for(std::vector<std::list<RangePair>>::iterator vect_it=vec_rp_list.begin(); vect_it!=vec_rp_list.end();vect_it++){
        vect_it->sort(RangePair::greaterScore);
        for (std::list<RangePair>::iterator i = vect_it->begin();
             i != vect_it->end(); i++)
            add_clean_included(i);
    }
}
//---------------------------------------------------------------------------
void BLRMatchAlign::write(std::ostream &out) {
    for (MapAlign::iterator m = map_align.begin(); m != map_align.end(); m++) {
        while (!m->second.empty()) {
            Key key = m->first;
            SDGString numQuery = SDGString(key.first);
            SDGString numSubject = SDGString(key.second);
            RangePair range_pair = m->second.back();
            SDGString numChr = SDGString(range_pair.getRangeQ().getNumChr());
            SDGString nameSeq = range_pair.getRangeQ().getNameSeq();
            m->second.pop_back();
            range_pair.write(out);
        }
    }

}