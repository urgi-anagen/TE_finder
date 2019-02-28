//
// Created by Hadi Quesneville on 2019-01-02.
//

#include "BLRMatchPath.h"
//----------------------------------------------------------------------------
void BLRMatchPath::insert(RangePairSet &rangePair)
//insert a RangePairSet in map_path at the right place
{
    std::list<RangePairSet> &al_list
            = map_path[rangePair.getRangeQ().getNumChr()];

    al_list.push_back(rangePair);
}
//---------------------------------------------------------------------------
void BLRMatchPath::read(const BLRJoinParameter& p, std::istream& input_path, int verbose)
{
    unsigned countseqS=0,countseqQ=0;

    name2numQ.clear();
    name2numS.clear();
    num2nameQ.clear();
    num2nameS.clear();
    para=p;

    unsigned long currentId=0, previousId=0;
    RangePairSet rps;

    while(input_path)
    {
        RangePairSet rp;
        rp.readtxt(input_path);
        if(input_path)
        {
            if(para.getEvalFilter()<rp.getE_value()
               || para.getIdFilter()>rp.getIdentity()
               || para.getLenFilter()>rp.getLength())
                continue;

            std::map<std::string,long>::iterator it
                    =name2numQ.find(rp.getRangeQ().getNameSeq());

            if(it==name2numQ.end()) //unknown name -> add to name list to associate a number
            {
                name2numQ[rp.getRangeQ().getNameSeq()]=++countseqQ;
                num2nameQ[countseqQ]=rp.getRangeQ().getNameSeq();
                rp.getRangeQ().setNumChr(countseqQ);
            }
            else
                rp.getRangeQ().setNumChr(it->second); // name already found -> add number


            it=name2numS.find(rp.getRangeS().getNameSeq());
            if(it==name2numS.end())
            {
                name2numS[rp.getRangeS().getNameSeq()]=++countseqS;
                num2nameS[countseqS]=rp.getRangeS().getNameSeq();
                rp.getRangeS().setNumChr(countseqS);
            }
            else
                rp.getRangeS().setNumChr(it->second);


            //remove sequence name to gain space
            rp.getRangeQ().setNameSeq("");
            rp.getRangeS().setNameSeq("");

            currentId = rp.getId();
            if (previousId != currentId && previousId!=0) //new path
            {
                //save previous path
                rps.setPath(para.getDist_pen(),0.0,para.getGap_pen());
                insert(rps);

                //clear for new path
                rps.clear();
            }
            rps.addPath(rp);
            previousId=currentId;
        }


    }
    rps.setPath(para.getDist_pen(),0.0,para.getGap_pen());
    insert(rps);

    if(verbose>0)
    {
        std::cout<<"nb of matches: "<<getNbMatchesInMapPath()<<std::endl;
        std::cout<<"nb of distinct queries: "<<getNbQseq()<<std::endl;
        std::cout<<"nb of distinct subjects: "<<getNbSseq()<<std::endl;
    }
}
//---------------------------------------------------------------------------
void BLRMatchPath::writeAttribute(std::ostream &out) {
    unsigned id=0;
    for (MapPath::iterator m = begin(); m != end(); m++) {
        m->second.sort(RangePair::less);
        for(std::list<RangePairSet>::iterator it=m->second.begin(); it != m->second.end(); it++){
            std::string query_name = num2nameQ[it->getNumQuery()];
            it->writeRpsAttr(out, ++id, query_name, num2nameS);
        }
    }
}
//---------------------------------------------------------------------------
void BLRMatchPath::write(std::ostream &out) {
    unsigned id=0;
    for (MapPath::iterator m = begin(); m != end(); m++) {
        m->second.sort(RangePair::less);
        for(std::list<RangePairSet>::iterator it=m->second.begin(); it != m->second.end(); it++){
            std::string query_name = num2nameQ[it->getNumQuery()];
            it->write(out, ++id, query_name, num2nameS);
        }
    }
}
//---------------------------------------------------------------------------
void BLRMatchPath::writeBED(const SDGString &filename, const SDGString &color) {
    std::ostringstream bedStream;
    writeBED(bedStream, color);
    std::ofstream bedFile(filename);
    bedFile << bedStream.str();
}

//---------------------------------------------------------------------------
void BLRMatchPath::writeBED(std::ostream &out, const SDGString &color) {

    for (MapPath::iterator m = begin(); m != end(); m++) {
        for(std::list<RangePairSet>::iterator it=m->second.begin(); it != m->second.end(); it++){
            std::string query_name = num2nameQ[it->getNumQuery()];
            it->writeBED(out, query_name, num2nameS, color);
        }
    }
}
//---------------------------------------------------------------------------
void BLRMatchPath::writeSeq(const SDGString &filename, int verbose) {
    if (verbose > 0)
        std::cout << "writing 'fasta' file..." << std::endl;

    SDGBioSeqDB query_db;
    if (para.getQuery() != "<not set>")
        query_db.load(para.getQuery());
    else{
        std::cerr<<"Error: Query databank not found!";
        exit(0);
    }

    SDGFastaOstream out(filename);
    int count_skip=0;

    for (MapPath::iterator iter_hash = begin();
         iter_hash != end(); iter_hash++) {

        std::string query_name = num2nameQ[iter_hash->first];
//        std::string subject_name = num2nameS[iter_hash->first.second];
//
//        std::string subseqname(subject_name, 0, subject_name.find(" "));
        SDGBioSeq chr=query_db.find(query_name);
        unsigned chr_len = chr.length();

        for (std::list<RangePairSet>::iterator i
                = iter_hash->second.begin(); i != iter_hash->second.end();
             i++) {
            if(i->getLength()<14)// temporaire !!??
            {
                count_skip++;
                continue;
            }
            long is=i->getRangeQ().getStart();
            long ie=i->getRangeQ().getEnd();
            if(is<=ie)
            {
                if(is>=chr_len)
                    continue;
            }
            else
            {
                if(ie>=chr_len)
                    continue;
            }
            SDGBioSeq s=newSDGMemBioSeq("");
            std::ostringstream name;
            if(i->getNbRangePairs()<=1)
            {
                if(is<=ie)
                {
                    s=chr.subseq(is-1,ie-is+1);
                }
                else
                {
                    s=chr.subseq(ie-1,is-ie+1);
                    s=s.complement();
                }
                name<<i->getRangeS().getNameSeq()<<" "<<i->getRangeQ().getNameSeq()
                    <<" {Fragment} "<<is<<".."<<ie;
            }
            else
            {
                SDGBioSeq sr=newSDGMemBioSeq("");
                name<<i->getRangeS().getNameSeq()<<" "<<i->getRangeQ().getNameSeq()
                    <<" {Fragment} ";
                bool first=true;
                for(std::list<RangePair>::const_iterator r= i->begin();
                    r!=i->end();r++)
                {
                    ulong rs=r->getRangeQ().getStart();
                    ulong re=r->getRangeQ().getEnd();
                    if(r->isPlusStrand())
                        sr=chr.subseq(rs-1,re-rs+1);
                    else
                    {
                        sr=chr.subseq(re-1,rs-re+1);
                        sr=sr.complement();
                    }
                    s+=sr;
                    if(!first) name<<",";
                    else first=false;
                    name<<rs<<".."<<re;
                }
            }
            s.setDE(name.str());
            out<<s;
        }
    }

    if(count_skip>0)
        std::cout<<count_skip<<" skipped!"<<std::endl;

    if (verbose > 0)
        std::cout << " done" << std::endl;
}
//---------------------------------------------------------------------------
void BLRMatchPath::writeMatch(const SDGString &filename, int verbose) {
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

    SDGBioSeqDB query_db,subject_db;
    if (para.getQuery() != "<not set>")
        query_db.load(para.getQuery());
    else{
        std::cerr<<"Error: Query databank not found!";
        exit(0);
    }
    if (para.getBank() != "<not set>")
        subject_db.load(para.getBank());
    else{
        std::cerr<<"Error: Subject databank not found!";
        exit(0);
    }

    unsigned path_id = 0;
    for (MapPath::iterator iter_hash = begin();
         iter_hash != end(); iter_hash++) {

        std::string query_name = num2nameQ[iter_hash->first];
        unsigned querylen = query_db.find(query_name).length();


        for (std::list<RangePairSet>::iterator iter_list
                = iter_hash->second.begin(); iter_list != iter_hash->second.end();
             iter_list++) {
            //iter_list->view();

            RangeAlignSet rasQ = iter_list->getRangeAlignSetQ();
            RangeAlignSet rasS = iter_list->getRangeAlignSetS();

            std::string subject_name = num2nameS[iter_list->getRangeS().getNumChr()];
            std::string subseqname(subject_name, 0, subject_name.find(" "));
            unsigned subjectlen = subject_db.find(subject_name).length();

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