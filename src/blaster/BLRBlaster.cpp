/**
 * \file BLRBlaster.cpp
 */

#include <sys/stat.h>
#include "BLRBlaster.h"

void BLRBlaster::run(int verbose) {
    if (verbose > 0) {
        std::cout << "subjects: " << para.getBankCut() << std::endl;
        std::cout << "queries: " << para.getQueryCut() << std::endl;
    }

    if (verbose > 0)
        std::cout << "loading query bank..." << std::flush;
    SDGBioSeqDB query_db(para.getQueryCut());
    if (verbose > 0)
        std::cout << " done" << std::endl;

    unsigned nbseq = para.getNumdo();
    for (unsigned i = 1; i <= nbseq; i++)
        waiting_seq.push_back(i);
    if (verbose > 0)
        std::cout << "nb of queries to process: " << waiting_seq.size() << std::endl;

    count_treated = 0;
    update_filename = para.getBlasterFileName() + ".seq_treated";
    std::ifstream file(update_filename);
    if (verbose > 0)
        std::cout << "queries already treated: ";
    while (file) {
        char buff[2048];
        unsigned num;
        file.getline(buff, 2048);
        if (!file)
            break;
        num = atoi(buff);
        waiting_seq.remove(num);
        count_treated++;
    }
    if (verbose > 0)
        std::cout << count_treated << std::endl;
    if (!count_treated)
        blasting->pressdb(verbose);

    if (para.getPrepare())
        return;

    if (waiting_seq.empty()) {
        std::cout << "All " << para.getType() << " are already done !" << std::endl;
        return;
    }

    std::list<unsigned> batch_num_seq;

    // for all remaining sequences
    while (!waiting_seq.empty()) {
        batch_num_seq.clear();

        //build batch of sequence to reach the total length of para.getLength()
        unsigned sum_len = 0;
        unsigned num_seq = 0, first_num_seq = 0;
        while ((para.getLength() == 0 || sum_len < para.getLength()) && !waiting_seq.empty()) {
            bool len_ok = false;
            SDGBioSeq s;
            while (!len_ok && !waiting_seq.empty()) {
                len_ok = true;
                num_seq = waiting_seq.front();
                if (sum_len == 0) first_num_seq = num_seq;
                s = query_db[num_seq - 1];
                if (s.length() < min_len) {
                    waiting_seq.pop_front();
                    len_ok = false;
                    if (verbose > 1)
                        std::cout << "skip sequence #" << num_seq
                                  << " < " << min_len << "bp" << std::endl;
                }
            }
            if (!waiting_seq.empty()) {
                batch_num_seq.push_back(num_seq);
                waiting_seq.pop_front();
                sum_len += s.length();
                if (verbose > 1) {
                    std::cout << "add sequence #" << num_seq
                              << " to batch #" << first_num_seq
                              << ", len=" << s.length();
                    std::cout << " sum=" << sum_len << std::endl;
                }
            }
        }
       if(waiting_seq.empty() && verbose > 1) {
           std::cout << "no more sequences !" << std::endl;
       }
        if (batch_num_seq.empty()){
            if (verbose > 1) {
                std::cout << "empty batch !" << std::endl;
            }
            break;
        }

        //run blast
        blasting->prepblast(batch_num_seq, query_db, first_num_seq);
        blasting->blast(verbose);

        //read results
        std::string matchfile = query_name + "_" + bank_name.afterlast("/")
                              + "_" + SDGString(first_num_seq) + ".res";
        BLRMatchList match_list;
        if (verbose > 1)
            std::cout << "read file '" << matchfile << "': ";
        match_list.read_blast_results(matchfile,
                                      para.get_is_wuBlast(),
                                      para.getEvalFilter(), para.getLenFilter(),para.getIdFilter() );
        if (verbose > 1)
            std::cout << match_list.getSize() << " HSPs" << std::endl;
        SDGString rm_cmd = "rm -f " + matchfile;
        system(rm_cmd);

        if (para.getBankCut() == para.getQueryCut()) //  remove_self_hits self hits when all-by-all
            match_list.remove_self_hits(para.getLength(),para.getOver());

        // save results in a file
        match_list.save_list(listfilenamecut, count_treated);
        count_treated++;
        if (verbose > 1)
            std::cout << "treated=" << count_treated << std::endl;
        update_seqtreatedfile(batch_num_seq); // update_seqtreatedfile seq_treated file list
    }

    if (verbose > 0)
        std::cout << "loading raw HSPs..." << std::endl;
    load_raw(verbose);
    if (verbose > 0)
        std::cout << "gluing HSPs..." << std::endl;
    glue(verbose);
    if (verbose > 0)
        std::cout << "clean self HSPs..." << std::endl;
    clean_self(verbose);
    if (verbose > 0)
        std::cout << "saving HSPs..." << std::endl;
    save_align(verbose);
    if (verbose > 0)
        std::cout << "Blaster was run." << std::endl;
    para.write(para.getParameterFileName());
}


void BLRBlaster::update_seqtreatedfile(const std::list<unsigned> &nseq) {
    std::ofstream file(update_filename, std::ios::app);
    if (file)
        for (std::list<unsigned>::const_iterator i = nseq.begin(); i != nseq.end(); i++)
            file << *i << std::endl;
    else
        std::cerr << "Can't update_seqtreatedfile batch " << nseq.front() << " !!" << std::endl;
}

void BLRBlaster::insert(RangePair &rangePair) {
    //insert rangePair in the right place
    std::list<RangePair> &al_list
            = map_align[Key(rangePair.getRangeQ().getNumChr(),
                            rangePair.getRangeS().getNumChr())];

    std::list<RangePair>::iterator r = al_list.begin();
    while (r != al_list.end() && *r < rangePair) r++;
    al_list.insert(r, rangePair);
}

void BLRBlaster::add_raw(const BlastMatch &align)
//add raw blast matches, compute absolute coordinate, no post processing
{
    if (para.getEvalFilter() < align.getE_value()
        || para.getIdFilter() > align.getIdentity()
        || para.getLenFilter() > align.getLength())
        return;

    RangePair rangePair(align);

    if (same_db) {
        if (rangePair.getRangeQ().getNumChr() > rangePair.getRangeS().getNumChr())
            return;
        if (rangePair.getRangeQ().getNumChr() == rangePair.getRangeS().getNumChr()
            && rangePair.getRangeQ().getMin() > rangePair.getRangeS().getMin())
            rangePair.invertQuerySubject();
    }

    long num_chr = rangePair.getRangeQ().getNumChr();
    RangeSeq range_seqdb = queryCut_db.getRange(num_chr);


    unsigned long start = range_seqdb.getStart();

    rangePair.getRangeQ().set(range_seqdb.getChr(),
                              name2numQ[range_seqdb.getChr()],
                              start + rangePair.getRangeQ().getStart() - 1,
                              start + rangePair.getRangeQ().getEnd() - 1);

    num_chr = rangePair.getRangeS().getNumChr();
    range_seqdb = subjectCut_db.getRange(num_chr);

    start = range_seqdb.getStart();

    if (same_db) {
        std::map<std::string, long>::iterator i = name2numQ.find(range_seqdb.getChr());
        if (i == name2numQ.end()) {
            std::cerr << "Chr: >>" << range_seqdb.getChr() << "<< not found!" << std::endl;
            throw SDGException(NULL, "Bank problem", -1);
        } else {
            rangePair.getRangeS().set(range_seqdb.getChr(), i->second,
                                      start + rangePair.getRangeS().getStart() - 1,
                                      start + rangePair.getRangeS().getEnd() - 1);
        }
    } else {
        std::map<std::string, long>::iterator i = name2numS.find(range_seqdb.getChr());
        if (i == name2numS.end()) {
            std::cerr << "Chr: >>" << range_seqdb.getChr() << "<< not found!" << std::endl;
            throw SDGException(NULL, "Bank problem", -1);
        } else {
            rangePair.getRangeS().set(range_seqdb.getChr(), i->second,
                                      start + rangePair.getRangeS().getStart() - 1,
                                      start + rangePair.getRangeS().getEnd() - 1);
        }
    }

    insert(rangePair);
}

void BLRBlaster::load_raw(int verbose)
// load from a raw blast matches file
{
    unsigned nbRaw = 0;
    queryCut_db.init(para.getQueryCut());
    subjectCut_db.init(para.getBankCut());

    same_db = (para.getBank() == para.getQuery());

    std::ifstream fileQ(para.getQuery());
    char buff[2048];

    unsigned long num = 0;
    while (fileQ) {
        fileQ.getline(buff, 2048);
        if (*(buff) == '>') {
            num++;
            SDGString chr_name(&buff[1]);
            chr_name = chr_name.trimL();
            chr_name = chr_name.trimR();
            name2numQ[chr_name] = num;
            num2nameQ[num] = chr_name;
        }
    }
    fileQ.close();

    if (!same_db) {
        std::ifstream fileS(para.getBank());
        num = 0;
        while (fileS) {
            fileS.getline(buff, 2048);
            if (*(buff) == '>') {
                num++;
                SDGString chr_name(&buff[1]);
                chr_name = chr_name.trimL();
                chr_name = chr_name.trimR();
                name2numS[chr_name] = num;
                num2nameS[num] = chr_name;
            }
        }
        fileS.close();
    }

    query_db.load(para.getQuery());
    subject_db.load(para.getBank());

    std::ifstream align_file(listfilenamecut);
    if (align_file.bad()) {
        std::cout << "BLRBlaster ERROR: " << listfilenamecut
                  << " could not be open" << std::endl;
        throw SDGException(NULL, "Blast results error", -1);
    }
    while (align_file) {
        BlastMatch align;
        align.readlst(align_file);
        if (align_file)
            add_raw(align);
        ++nbRaw;
    }
    align_file.close();
    if (verbose > 0)
        std::cout << "loaded " << nbRaw << " HSPs" << std::endl;
}

void BLRBlaster::save_align(int verbose) {
    unsigned nbSaved = 0;
    std::ofstream align_file(para.getParameterFileName().beforelast(".param") + ".align");
    if (align_file.bad()) {
        std::cout << "BLRBlaster ERROR: " << para.getParameterFileName().beforelast(".param") << " could not be open"
                  << std::endl;
        throw SDGException(NULL, "fatal error", -1);
    }
    for (MapAlign::iterator m = map_align.begin(); m != map_align.end(); m++) {
        for (std::list<RangePair>::iterator r = m->second.begin(); r != m->second.end(); r++) {
            if (para.getAll_by_all())
                if (r->getRangeQ().getNameSeq() == r->getRangeS().getNameSeq()
                    && r->getRangeQ().getStart() == r->getRangeS().getStart()
                    && r->getRangeQ().getEnd() == r->getRangeS().getEnd())
                    continue;
            //r->view();
            r->write(align_file);
            ++nbSaved;
        }
    }
    align_file.close();
    if (verbose > 0)
        std::cout << "saved " << nbSaved << " HSPs" << std::endl;
}

void BLRBlaster::view_align(void) {
    for (MapAlign::iterator m = map_align.begin(); m != map_align.end(); m++) {
        for (std::list<RangePair>::iterator r = m->second.begin(); r != m->second.end();
             r++) {
            r->view();
        }
    }
}

void BLRBlaster::clean_self(int verbose)
// remove_self_hits matches when query fragments are the same
{
    if (para.getBank() != para.getQuery()) return;
    for (MapAlign::iterator m = map_align.begin(); m != map_align.end(); m++) {
        if (m->first.first == m->first.second) {
            std::list<RangePair> &al_list = m->second;
            for (std::list<RangePair>::iterator r = al_list.begin();
                 r != al_list.end(); r++) {
                if (r->getRangeQ().getStart() == r->getRangeS().getStart()
                    && r->getRangeQ().getEnd() == r->getRangeS().getEnd()) {
                    r = al_list.erase(r);
                    r--;
                }
            }
        }
    }
}

void BLRBlaster::glue(int verbose)
//  glue matches stopped by cutting of bank
{
    for (MapAlign::iterator m = map_align.begin(); m != map_align.end(); m++) {
        std::list<RangePair> &al_list = m->second;
        al_list.sort(RangePair::less);
        std::list<RangePair>::iterator r2, r1 = al_list.begin();
        while (r1 != al_list.end()) {
            r2 = r1;
            r2++;
            while (r2 != al_list.end() && r1->overlap(*r2) && r1 != r2) {
                if (r1->getRangeQ().isPlusStrand()
                    == r2->getRangeQ().isPlusStrand()
                    && r1->getRangeS().isPlusStrand()
                       == r2->getRangeS().isPlusStrand()
                        ) {
                    unsigned overq, overs;
                    if (r1->getRangeQ() < r2->getRangeQ())
                        overq = r1->getRangeQ().getMax() - r2->getRangeQ().getMin();
                    else
                        overq = r2->getRangeQ().getMax() - r1->getRangeQ().getMin();
                    if (r1->getRangeS() < r2->getRangeS())
                        overs = r1->getRangeS().getMax() - r2->getRangeS().getMin();
                    else
                        overs = r2->getRangeS().getMax() - r1->getRangeS().getMin();
                    float ratio = (float) overq / overs;
                    if (ratio > 0.95 && ratio < 1.05) {
                        r1->merge(*r2);
                        r2 = al_list.erase(r2);
                        r2--;
                    }
                } //test strand
                r2++;
            } //overlap loop
            r1++;
        } //r1 loop
    }
}

void BLRBlaster::cleanTmpFiles(int verbose) {
    std::ostringstream cmd;
    struct stat stFileInfo;
    int intStat;

    std::list<std::string> lFiles;
    lFiles.push_back("formatdb.log");
    lFiles.push_back("last_time_stamp.log");
    lFiles.push_back("error.log");
    lFiles.push_back(para.getBlasterFileName() + ".raw");
    lFiles.push_back(para.getBlasterFileName() + ".seq_treated");
    lFiles.push_back(para.getQuery() + ".Nstretch.map");
    lFiles.push_back(para.getQueryCut());
    lFiles.push_back(para.getBank() + ".Nstretch.map");
    lFiles.push_back(para.getBankCut());

    if (!para.get_is_wuBlast()) {
        if (para.getType() == "blastn" || para.getType() == "megablast") {
            lFiles.push_back(para.getBankCut() + ".nhr");
            lFiles.push_back(para.getBankCut() + ".nin");
            lFiles.push_back(para.getBankCut() + ".nsq");
        } else {
            lFiles.push_back(para.getBankCut() + ".phr");
            lFiles.push_back(para.getBankCut() + ".pin");
            lFiles.push_back(para.getBankCut() + ".psq");
        }
    }

    if (para.get_is_wuBlast()) {
        if (para.getType() == "blastn") {
            lFiles.push_back(para.getBankCut() + ".xnd");
            lFiles.push_back(para.getBankCut() + ".xns");
            lFiles.push_back(para.getBankCut() + ".xnt");
        } else {
            lFiles.push_back(para.getBankCut() + ".xpd");
            lFiles.push_back(para.getBankCut() + ".xps");
            lFiles.push_back(para.getBankCut() + ".xpt");
        }
    }

    for (std::list<std::string>::iterator itFile = lFiles.begin(); itFile != lFiles.end(); ++itFile) {
        intStat = stat((*itFile).c_str(), &stFileInfo);
        if (intStat == 0) {
            cmd << "rm " << (*itFile).c_str();
            system(cmd.str().c_str());
            cmd.str("");
        }
    }
}
