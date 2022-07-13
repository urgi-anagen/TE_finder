#include <SDGBioSeqDB.h>
#include <FragJoin.h>
#include <Lalign.h>
#include <FastLalign.h>
#include <FastExtAlign.h>
#include "Hasher.h"

//-------------------------------------------------------------------------
// Search for diagonal of word matches with distance
void Hasher::diagSearchDist(unsigned numseqQ, Diag_map &diag_map,
                            unsigned connect_dist, unsigned kmer_size, unsigned min_frag_size,
                            std::list< RangePair >& frag, unsigned verbose) {

    unsigned count_frag = 0;
    unsigned curr_seq = 0;

    for (auto &iter_seq : diag_map) { // iter seq
        if (iter_seq.size() > 2) {
            bool extending = false;
            iter_seq.sort();
            auto iter_diag = iter_seq.begin();
            Diag prev_d = *iter_diag;
            unsigned start = 0;
            unsigned end = 0;
            int diag = 0;
            unsigned score = 0;

            while (++iter_diag != iter_seq.end()) {
                Diag curr_d = *iter_diag;
                curr_seq = curr_d.wpos.numSeq;

                if (prev_d.diag == curr_d.diag
                    && prev_d.wpos.numSeq == curr_d.wpos.numSeq
                    && (prev_d.wpos.pos + connect_dist >= curr_d.wpos.pos
                        && prev_d.wpos.pos + diag + connect_dist >= curr_d.wpos.pos + diag)
                        ) {
                    if (extending) //extending
                    {
                        end = curr_d.wpos.pos;
                        score++;
                    } else //first hit (2 kmers found at correct distance)
                    {
                        diag = prev_d.diag;
                        start = prev_d.wpos.pos;
                        end = curr_d.wpos.pos;
                        extending = true;
                        score = 1;
                    }
                } else //stop extension if distance between kmer too long
                if (extending) {
                    if (end + kmer_size - start - 1 >= min_frag_size) {
                        count_frag++;
                        frag.push_back(record_frag(start, end, diag,
                                                   score, numseqQ, curr_seq, count_frag));
                    }
                    extending = false;
                }
                prev_d = curr_d;
            } //end while
            if (extending) // Record hit at the end of the loop
            {
                if (end + kmer_size - start - 1 >= min_frag_size) {
                    count_frag++;
                    frag.push_back(record_frag(start, end, diag,
                                               score, numseqQ, curr_seq, count_frag));
                }
            }
        } //end size>2,
    } // seq loop

    if (verbose > 0) {
        std::cout << "Fragments number founds:" << count_frag << std::endl;
    }
}
//-------------------------------------------------------------------------
// Search for Alignments
void Hasher::search(const BioSeq& sequence, unsigned start, unsigned end, unsigned numseq, unsigned connect_dist,
		unsigned min_frag_size, bool repeat, std::list< RangePair >& frag, unsigned verbose)
{
    if(end-start+1<min_frag_size) return;
	clock_t clock_begin, clock_end;
	clock_begin = clock();
	std::cout<<"hashing query sequence #"<<numseq<<" from "<<start<<" to "<<end<<" ..."<<std::flush;

    Diag_map diag_map(subject_names.size()+1);
    if(algorithm==0)
        matchKmers(sequence, start, end, repeat, diag_map);
    else if(algorithm==1)
        matchKmersHole(sequence, start, end, repeat, diag_map);
    else if(algorithm==2)
        matchKmersMinimizer(sequence, start, end, repeat, diag_map);

	clock_end = clock();
	std::cout<<" --> Time spent: "<<(double)(clock_end-clock_begin)/CLOCKS_PER_SEC<<" seconds"<<std::endl;

	clock_begin = clock();
	std::cout<<"run_test_search_wSW fragments..."<<std::flush;
    diagSearchDist(numseq, diag_map, connect_dist, kmer_size, min_frag_size, frag, verbose-1);

	diag_map.clear();

	std::cout<<"ok"<<std::endl;
	clock_end = clock();
	std::cout<<" --> Time spent: "<<(double)(clock_end-clock_begin)/CLOCKS_PER_SEC<<" seconds"<<std::endl;
}
//-------------------------------------------------------------------------
// merge found fragments
void Hasher::fragJoin(std::list< RangePair >& frag) const
{
    //separate fragments on same query and same subject sequences
    std::map<std::pair<unsigned,unsigned>, std::list<RangePair> > map_frag;
    auto it=frag.begin();
    while(it!=frag.end())
    {
        std::pair<unsigned,unsigned> key(it->getRangeQ().getNumChr(),it->getRangeS().getNumChr());
        map_frag[key].push_back(*it);
        it=frag.erase(it);
    }
    double gap_pen=pen_join;
    double dist_pen=pen_join;
    FragJoin fragJoin(dist_pen, 0, gap_pen,0);

   // merge contiguous fragments
    for (auto & item : map_frag) {
        std::list<RangePairSet> jfrag;
        fragJoin.align_all(item.second, jfrag);
        for(auto & f : jfrag){
            frag.emplace_back(f);
        }
    }
}
//-------------------------------------------------------------------------
// Merge overlapping rangePair
void Hasher::fragMerge(const std::list< RangePair >& frag, std::list< RangePair >& frag_merge)
{
    frag_merge=frag;
    frag_merge.sort(RangePair::less);

    unsigned size=frag_merge.size();
    if(size>=2){
        auto curr_frag_it=frag_merge.begin();
        auto next_frag_it =curr_frag_it;
        next_frag_it++;

        while(next_frag_it != frag_merge.end()) {
            if (curr_frag_it->overlapQ(*next_frag_it)) {
                //TODO chose a subject name if different
                curr_frag_it->merge(*next_frag_it);
                next_frag_it = frag_merge.erase(next_frag_it);
            } else{
                curr_frag_it++;
                next_frag_it++;
            }
        }
    }
}
//-------------------------------------------------------------------------
// Merge overlapping rangePair
unsigned Hasher::fragCoverage(const std::list< RangePair >& frag)
{
    unsigned coverage=0;
    std::list< RangePair > frag_merge;
    fragMerge(frag,frag_merge);

    for(auto & it : frag_merge)
        coverage+=it.getLength();
    return coverage;
}
//-------------------------------------------------------------------------
// Stats on rangePair lists
unsigned Hasher::fragScoreIdentityStat(const std::list< RangePair >& frag, double quantile, unsigned& coverage)
{
    coverage=0;
    if(frag.empty()){
        return 0;
    }
    std::vector<unsigned> score_list;
    std::vector<double> identity_list;
    for(const auto & curr_frag_it : frag) {
        unsigned len=curr_frag_it.getLength();
        coverage+=len;
        score_list.push_back(curr_frag_it.getScore());
        identity_list.push_back(curr_frag_it.getIdentity());
    }
    sort(score_list.begin(), score_list.end());
    unsigned nb_frag=score_list.size();
    unsigned min_score=score_list.front();
    unsigned max_score=score_list.back();
    unsigned qval_score=score_list[(int)std::floor((double)score_list.size() * quantile)];
    std::cout << "Frag number=" << nb_frag << " / "
              << "min score=" << min_score << " / "
              << "max score=" << max_score << " / "
              << "quantile score (" << quantile << ")=" << qval_score
             <<std::endl;
    sort(identity_list.begin(), identity_list.end());
    unsigned min_identity=identity_list.front();
    unsigned max_identity=identity_list.back();
    unsigned qval_identity=identity_list[(int)std::floor((double)identity_list.size() * quantile)];
    std::cout << "     min identity=" << min_identity << " / "
              << "max identity=" << max_identity << " / "
              << "quantile identity (" << quantile << ")=" << qval_identity
              <<std::endl;
    return qval_score;
}
//-------------------------------------------------------------------------
unsigned Hasher::fragLengthStat(const std::list< RangePair >& frag, double quantile)
{
    if(frag.empty()){
        return 0;
    }
    std::vector<unsigned> length_list;
    for(const auto & curr_frag_it : frag) {
        length_list.push_back(curr_frag_it.getLength());
    }
    sort(length_list.begin(), length_list.end());
    unsigned nb_frag=length_list.size();
    unsigned min_score=length_list.front();
    unsigned max_score=length_list.back();
    unsigned qval=length_list[(int)std::floor((double)length_list.size() * quantile)];
    std::cout << "Frag number=" << nb_frag << " / "
              << "min length=" << min_score << " / "
              << "max length=" << max_score << " / "
              <<"quantile length ("<<quantile<<")="<<qval
              <<std::endl;
    return qval;
}
//-------------------------------------------------------------------------
// Filter length on rangePair lists
void Hasher::fragLenFilter(std::list< RangePair >& frag, unsigned min_len)
{
    std::cout<<"--Filter fragments length <"<<min_len<<" ... "<<std::flush;
    auto frag_it=frag.begin();
    while (frag_it != frag.end()) {
        if (frag_it->getLength() < min_len) {
            frag_it = frag.erase(frag_it);
        } else { frag_it++; }
    }
    std::cout<<"done !"<<std::endl;
}
//-------------------------------------------------------------------------
// Filter score on rangePair lists
void Hasher::fragScoreIdentityFilter(std::list<RangePair> &frag, unsigned min_score, unsigned min_identity)
{
    std::cout<<"--Filter fragments score <"<<min_score<<" and identity <"<<min_identity<<" ... "<<std::flush;
    auto frag_it=frag.begin();
    while(frag_it != frag.end()) {
        if(frag_it->getScore()<min_score || frag_it->getIdentity()<min_identity){
            frag_it = frag.erase(frag_it);
        }else{frag_it++;}
    }
    std::cout<<"done !"<<std::endl;
}
//-------------------------------------------------------------------------
// Set rangePair score and identity
void Hasher::fragSeqAlign(std::list< RangePair >& frag,
                          const SDGString& fasta_queryfilename, const SDGString& fasta_subjectfilename,
                          bool reverse, unsigned verbose)
{
    FastaIstream query_in(fasta_queryfilename);
    if (!query_in) {
        std::cerr << "file:" << fasta_queryfilename << " does not exist!" << std::endl;
    }

    std::vector<BioSeq> subject_db;
    FastaIstream subject_in(fasta_subjectfilename);
    if (!subject_in) {
        std::cerr << "file:" << fasta_queryfilename << " does not exist!" << std::endl;
    }
    while (subject_in) {
        BioSeq seq;
        if (subject_in) {
            subject_in >> seq;
            subject_db.push_back(seq);
        }
    }

    unsigned numseq=0;
    while (query_in) {
        BioSeq seq;
        if (query_in)
            query_in >> seq;
        if(reverse){
            seq=seq.reverse();
        }
        numseq++;
        if(verbose>0) std::cout << seq.header << " len:" << seq.size() << " read!" << std::endl;
        for (auto & curr_frag_it : frag) {
            if (curr_frag_it.getRangeQ().getNumChr() == numseq) {
                if(verbose>0) curr_frag_it.view();
                // RangePair on the current query sequence
                BioSeq qseq;
                qseq = seq.subseq(curr_frag_it.getRangeQ().getMin()-1,
                                  curr_frag_it.getRangeQ().getLength());
                if (!curr_frag_it.getRangeQ().isPlusStrand()) {
                    qseq = qseq.complement();
                }
                unsigned qlen=qseq.size();

                if(verbose>0) std::cout << "query:  " << qseq << "-" << qlen << std::endl;
                // Get subject sequence

                BioSeq sseq = subject_db[curr_frag_it.getRangeS().getNumChr()-1];

                BioSeq fragsseq;
                fragsseq = sseq.subseq(curr_frag_it.getRangeS().getMin()-1,
                                       curr_frag_it.getRangeS().getLength());
                if (!curr_frag_it.getRangeS().isPlusStrand()) {
                    fragsseq = fragsseq.complement();
                }
                unsigned slen=fragsseq.size();
                if(verbose>0) std::cout << "subject:" << fragsseq << "-" << slen << std::endl;
                unsigned count=0;

                unsigned len_align=std::min(qlen,slen);
                for(unsigned i=0;i<len_align;i++)
                {
                    if(qseq[i]==fragsseq[i]) count++;
                }

                curr_frag_it.setIdentity(((double)count)/len_align*100);
                curr_frag_it.setScore(count);
                if(verbose>0) std::cout << "Score = " << count<<" identity = " << ((double)count)/len_align * 100<< std::endl;
            }
        }
    }
}
//-------------------------------------------------------------------------
// SW Align rangePair
void Hasher::fragSeqSWAlign(std::list< RangePair >& frag, unsigned ext_len,
                            const SDGString& fasta_queryfilename, const SDGString& fasta_subjectfilename,
                            bool reverse, unsigned verbose)
{
    //Lalign lalign(1);
    FastLalign lalign;
    lalign.setMismatch(5,4);

    FastaIstream query_in(fasta_queryfilename);
    if (!query_in) {
        std::cerr << "file:" << fasta_queryfilename << " does not exist!" << std::endl;
    }

    std::vector<BioSeq> subject_db;
    FastaIstream subject_in(fasta_subjectfilename);
    if (!subject_in) {
        std::cerr << "file:" << fasta_queryfilename << " does not exist!" << std::endl;
    }
    while (subject_in) {
        BioSeq seq;
        if (subject_in) {
            subject_in >> seq;
            subject_db.push_back(seq);
        }
    }

    unsigned numseq=0;
    unsigned nb_align=0;
    while (query_in) {
        BioSeq qseq;
        if (query_in)
            query_in >> qseq;
        if(reverse){
            qseq=qseq.reverse();
        }
        numseq++;
        if(verbose>0) std::cout << qseq.header << " len:" << qseq.size() << " read!" << std::endl;
        for (auto & curr_frag_it : frag) {
            if (curr_frag_it.getRangeQ().getNumChr() == numseq) {
                if(verbose>0) {
                    curr_frag_it.view(); }
                // RangePair on the current query sequence
                BioSeq fragqseq;
                unsigned qstart,qlen;

                if(curr_frag_it.getRangeQ().getMin()-1<ext_len)
                    qstart=0;
                else
                    qstart=curr_frag_it.getRangeQ().getMin()-1-ext_len;

                if(curr_frag_it.getRangeQ().getLength()+2*ext_len > qseq.size()){
                    qlen= qseq.size() - qstart;
                }
                else{
                    qlen=curr_frag_it.getRangeQ().getLength()+2*ext_len;
                }

                fragqseq = qseq.subseq(qstart, qlen);

                if (!curr_frag_it.getRangeQ().isPlusStrand()) {
                    fragqseq = fragqseq.complement();
                }
                qlen=fragqseq.size();

                // Get subject sequence

                BioSeq sseq = subject_db[curr_frag_it.getRangeS().getNumChr()-1];

                BioSeq fragsseq;
                unsigned sstart,slen;

                if(curr_frag_it.getRangeS().getMin()-1<ext_len)
                    sstart=0;
                else
                    sstart=curr_frag_it.getRangeS().getMin()-1-ext_len;

                if(curr_frag_it.getRangeS().getLength()+2*ext_len > sseq.size()){
                    slen=sseq.size()-sstart;
                }
                else{
                    slen=curr_frag_it.getRangeS().getLength()+2*ext_len;
                }


                fragsseq = sseq.subseq(sstart,slen);
                if (!curr_frag_it.getRangeS().isPlusStrand()) {
                    fragsseq = fragsseq.complement();
                }

                slen=fragsseq.size();

                if(verbose>0) std::cout<<"aligning sequence of length : "<<qlen<<" ... "<<std::flush;
                if(++nb_align % 100 ==0) std::cout<<std::endl;
                std::cout<<"."<<std::flush;
                lalign.setSeq(fragqseq, fragsseq);
                lalign.align();
                if(verbose>0) std::cout<<" done !"<<std::endl;
                if(verbose>1) lalign.view();

                if(lalign.getScore()>0){
                    curr_frag_it.setIdentity(lalign.getIdentity()*100);
                    curr_frag_it.setScore(lalign.getScore());
                    curr_frag_it.setLength(lalign.getLength());

                    if(curr_frag_it.getRangeQ().isPlusStrand()){
                        curr_frag_it.getRangeQ().setStart(lalign.getStartSeq1()+qstart);
                        curr_frag_it.getRangeQ().setEnd(lalign.getEndSeq1()+qstart);
                    }else{
                        curr_frag_it.getRangeQ().setEnd(fragqseq.size() - lalign.getEndSeq1() + qstart + 1);
                        curr_frag_it.getRangeQ().setStart(fragqseq.size() - lalign.getStartSeq1() + qstart + 1);
                    }
                    curr_frag_it.getRangeS().setStart(lalign.getStartSeq2()+sstart);
                    curr_frag_it.getRangeS().setEnd(lalign.getEndSeq2()+sstart);
                }else{
                    std::cout << "No SW alignment found !"<< std::endl;
                    if(verbose>0) {
                        std::cout << "query sequence ("<<qstart<<","<<qlen<<" len=" << qseq.size() << ") =" << fragqseq << std::endl;
                        std::cout << "subject sequence ("<<sstart<<","<<slen<<" len="<<sseq.size()<< ") =" << fragsseq << std::endl;
                    }
                }

                if(verbose>0) {
                    std::cout << "----result---->";
                    curr_frag_it.view(); }
                if(verbose>0) std::cout << "Score = " << curr_frag_it.getScore()<<" identity = " << curr_frag_it.getIdentity()<< std::endl;
            }
        }
    }
}
//-------------------------------------------------------------------------
// Extend boundaries by alignment
void Hasher::fragSeqExtAlign(std::list< RangePair >& frag, unsigned ext_len,
                            const SDGString& fasta_queryfilename, const SDGString& fasta_subjectfilename,
                            bool reverse, unsigned verbose)
{
    FastExtAlign align;
    align.setMismatch(5,4);

    FastaIstream query_in(fasta_queryfilename);
    if (!query_in) {
        std::cerr << "file:" << fasta_queryfilename << " does not exist!" << std::endl;
    }

    std::vector<BioSeq> subject_db;
    FastaIstream subject_in(fasta_subjectfilename);
    if (!subject_in) {
        std::cerr << "file:" << fasta_queryfilename << " does not exist!" << std::endl;
    }
    while (subject_in) {
        BioSeq seq;
        if (subject_in) {
            subject_in >> seq;
            subject_db.push_back(seq);
        }
    }

    unsigned numseq=0;
    unsigned nb_align=0;
    while (query_in) {
        BioSeq qseq;
        if (query_in)
            query_in >> qseq;
        if(reverse){
            qseq=qseq.reverse();
        }
        numseq++;
        if(verbose>0) std::cout << qseq.header << " len:" << qseq.size() << " read!" << std::endl;
        for (auto & curr_frag_it : frag) {
            if (curr_frag_it.getRangeQ().getNumChr() == numseq) {
                if(verbose>0) {
                    curr_frag_it.view(); }
                // RangePair on the current query sequence
                BioSeq fragqseq;
                unsigned qstart,qlen;

                if(curr_frag_it.getRangeQ().getMin()-1<ext_len)
                    qstart=0;
                else
                    qstart=curr_frag_it.getRangeQ().getMin()-1-ext_len;

                if(curr_frag_it.getRangeQ().getLength()+2*ext_len > qseq.size()){
                    qlen= qseq.size() - qstart;
                }
                else{
                    qlen=curr_frag_it.getRangeQ().getLength()+2*ext_len;
                }

                fragqseq = qseq.subseq(qstart, qlen);

                if (!curr_frag_it.getRangeQ().isPlusStrand()) {
                    fragqseq = fragqseq.complement();
                }
                qlen=fragqseq.size();

                // Get subject sequence

                BioSeq sseq = subject_db[curr_frag_it.getRangeS().getNumChr()-1];

                BioSeq fragsseq;
                unsigned sstart,slen;

                if(curr_frag_it.getRangeS().getMin()-1<ext_len)
                    sstart=0;
                else
                    sstart=curr_frag_it.getRangeS().getMin()-1-ext_len;

                if(curr_frag_it.getRangeS().getLength()+2*ext_len > sseq.size()){
                    slen=sseq.size()-sstart;
                }
                else{
                    slen=curr_frag_it.getRangeS().getLength()+2*ext_len;
                }


                fragsseq = sseq.subseq(sstart,slen);
                if (!curr_frag_it.getRangeS().isPlusStrand()) {
                    fragsseq = fragsseq.complement();
                }

                slen=fragsseq.size();

                if(verbose>0) std::cout<<"aligning sequence of length : "<<qlen<<" ... "<<std::flush;
                if(++nb_align % 100 ==0) std::cout<<std::endl;
                std::cout<<"."<<std::flush;
                align.setSeq(fragqseq, fragsseq);
                align.setStart(fragqseq.size()-ext_len,fragsseq.size()-ext_len,ext_len);
                unsigned score_dir=align.extend_dir(curr_frag_it.getScore());
                unsigned end_seq1=align.getEndSeq1();
                unsigned end_seq2=align.getEndSeq2();

                align.reset_align();
                align.setSeq(fragqseq, fragsseq);

                align.setStart(ext_len,ext_len,ext_len);
                unsigned score_rev=align.extend_rev(curr_frag_it.getScore());
                unsigned start_seq1=align.getEndSeq1();
                unsigned start_seq2=align.getEndSeq2();
                if(verbose>0) std::cout<<" done !"<<std::endl;

                curr_frag_it.setScore(curr_frag_it.getScore()+score_dir+score_rev);

                if(curr_frag_it.getRangeQ().isPlusStrand()){
                    curr_frag_it.getRangeQ().setStart(start_seq1+qstart);
                    curr_frag_it.getRangeQ().setEnd(end_seq1+qstart);
                }else{
                    curr_frag_it.getRangeQ().setEnd(fragqseq.size() - end_seq1 + qstart + 1);
                    curr_frag_it.getRangeQ().setStart(fragqseq.size() - start_seq1 + qstart + 1);
                }
                curr_frag_it.getRangeS().setStart(start_seq2+sstart);
                curr_frag_it.getRangeS().setEnd(end_seq2+sstart);


                if(verbose>0) {
                    std::cout << "----result---->";
                    curr_frag_it.view(); }
                if(verbose>0) std::cout << "Score = " << curr_frag_it.getScore()<<" identity = " << curr_frag_it.getIdentity()<< std::endl;
            }
        }
    }
}
//-------------------------------------------------------------------------
// Write a rangePair lists
void Hasher::fragAlignWrite(std::list< RangePair >& frag, const SDGString& qfilename, const SDGString& sfilename, std::ostream& out)
{
    std::vector<std::string> num2nameQ,num2nameS;

    std::ifstream fileQ(qfilename);
    std::string line;
    if (fileQ.is_open()) {
        while (std::getline(fileQ, line)) {
            if (line[0] == '>') {
                line.erase(0, 1);
                line.erase(line.find_last_not_of("\t\n\v\f\r ") + 1);
                num2nameQ.push_back(line);
            }
        }
    }
    fileQ.close();

    std::ifstream fileS(sfilename);
    if (fileS.is_open()) {
        while (std::getline(fileS, line)) {
            if (line[0] == '>') {
                line.erase(0, 1);
                line.erase(line.find_last_not_of("\t\n\v\f\r ") + 1);
                num2nameS.push_back(line);
            }
        }
    }
    fileS.close();

    for(auto & curr_frag_it : frag) {
        if(curr_frag_it.getNumQuery()>(long)num2nameQ.size()){
            std::cerr<<"Error query sequence "<<curr_frag_it.getRangeQ().getNameSeq()<<" number "<<curr_frag_it.getNumQuery()<<" doesn't exist!"<<std::endl;
            exit(EXIT_FAILURE);
        }
        SDGString qname=num2nameQ[curr_frag_it.getNumQuery()-1];

        if(curr_frag_it.getNumSubject()>(long)num2nameS.size()){
            std::cerr<<"Error subject sequence "<<curr_frag_it.getRangeS().getNameSeq()<<" number "<<curr_frag_it.getNumSubject()<<" doesn't exist!"<<std::endl;
            exit(EXIT_FAILURE);
        }
        SDGString sname=num2nameS[curr_frag_it.getNumSubject()-1];
        curr_frag_it.setQSName(qname,sname);
        curr_frag_it.write(out);
    }
}
//-------------------------------------------------------------------------
// Write a rangePair sequences
void Hasher::fragSeqWrite(const std::list< RangePair >& frag, const SDGString& fasta_filename, FastaOstream& out)
{
    FastaIstream in(fasta_filename);
    if (!in) {
        std::cerr << "file:" << fasta_filename << " does not exist!" << std::endl;
    }
    unsigned numseq=0;
    while (in) {
        BioSeq seq;
        if (in)
            in >> seq;
        numseq++;
        std::cout << seq.header << " len:" << seq.size() << " read!" << std::endl;
        for (const auto & curr_frag : frag) {
            if (curr_frag.getRangeQ().getNumChr() == numseq) {
                BioSeq sseq;
                if (curr_frag.getRangeQ().isPlusStrand()) {
                    sseq = seq.subseq(curr_frag.getRangeQ().getStart(), curr_frag.getRangeQ().getEnd() - curr_frag.getRangeQ().getStart() + 1);
                } else {
                    sseq = seq.subseq(curr_frag.getRangeQ().getEnd(), curr_frag.getRangeQ().getStart() - curr_frag.getRangeQ().getEnd() + 1);
                    sseq = sseq.complement();
                }
                std::istringstream subject_name(curr_frag.getRangeS().getNameSeq());
                std::string prefix_name;
                subject_name>>prefix_name;
                std::ostringstream name;
                name << prefix_name << " " << seq.header << ":"
                     << curr_frag.getRangeQ().getStart() << ".."
                     << curr_frag.getRangeQ().getEnd();
                sseq.header=name.str();
                out << sseq;
            }
        }
    }
}
//-------------------------------------------------------------------------
// Write a rangePair sequences
void Hasher::fragMergeSeqWrite(const std::list< RangePair >& frag, const SDGString& fasta_filename, FastaOstream& out)
{
    std::list< RangePair > frag_merge;
    fragMerge(frag,frag_merge);
    fragSeqWrite(frag_merge,fasta_filename,out);
}