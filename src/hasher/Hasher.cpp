#include <SDGBioSeqDB.h>
#include <FragJoin.h>
#include "Hasher.h"

//-------------------------------------------------------------------------
// Search for diagonal of word matches with distance
void Hasher::diagSearchDist(unsigned numseqQ, std::vector<std::list<Diag> > &diag_map,
                            unsigned connect_dist, unsigned kmer_size, unsigned min_frag_size,
                            std::list< RangePair >& frag, unsigned verbose) {

    unsigned count = 0;
    unsigned curr_seq = 0;
    for (auto &iter_seq : diag_map) {
        count = 0;
        unsigned size = iter_seq.size();
        if (size > 2) {
            iter_seq.sort();
            bool extending=false;
            unsigned start = 0;
            unsigned end = 0;
            unsigned score = 0;
            long diag = 0;

            auto iter_diag = iter_seq.begin();
            Diag prev_d = *iter_diag;
            while (++iter_diag != iter_seq.end()) {
                Diag curr_d = *iter_diag;
                curr_seq = curr_d.wpos.numSeq;

                if (prev_d.diag == curr_d.diag
                    && prev_d.wpos.numSeq == curr_d.wpos.numSeq
                    && prev_d.wpos.pos + connect_dist >= curr_d.wpos.pos) {
                    if (extending) //extending
                    {
                        end = curr_d.wpos.pos;
                        score++;
                    } else //first hit (2 kmers found at correct distance)
                    {
                        diag = prev_d.diag;
                        start = prev_d.wpos.pos;
                        end = curr_d.wpos.pos;
                        extending=true;
                        score = 1;
                    }
                } else //stop extension if distance between kmer too long
                if (extending) {
                    if (end + kmer_size - start - 1 >= min_frag_size) {
                        count++;
                        frag.push_back(record_frag(start, end, diag,
                                                   score, numseqQ, curr_seq, count));
                    }
                    extending=false;
                }
                prev_d = curr_d;
            } //end for
            if (extending) // Record hit at the end of the loop
            {
                if (end + kmer_size - start - 1 >= min_frag_size) {
                    count++;
                    frag.push_back(record_frag(start, end, diag,
                                               score, numseqQ, curr_seq, count));
                }
            }
        } //end size>2, diag_map loop
    }

    if (verbose > 0) {
        std::cout << "Fragments number founds:" << count << std::endl;
    }
}
//-------------------------------------------------------------------------
// Search for diagonal of word matches with score and penalty
void Hasher::diagSearchScore(unsigned numseqQ, std::vector<std::list<Diag> > &diag_map,
                             unsigned min_frag_size,
                             std::list<RangePair> &frag, unsigned verbose) {
    double penalty = (double)1/wdist;
    unsigned count;

    for (auto &iter_seq : diag_map) { // iter diagonals
        unsigned size = iter_seq.size();
        count=0;
        if (size > 2) {
            iter_seq.sort();
            bool extending = false;
            unsigned start = 0;
            unsigned end = 0;
            unsigned score = kmer_size;
            long diag = 0;

            auto iter_diag = iter_seq.begin();
            Diag curr_d, prev_d = *iter_diag;
            while (++iter_diag != iter_seq.end()) {
                curr_d = *iter_diag;
                if (prev_d.diag == curr_d.diag && prev_d.wpos.numSeq == curr_d.wpos.numSeq) {
                    long extended_score;
                    long dist = curr_d.wpos.pos - prev_d.wpos.pos - kmer_size;
                    if (dist < 0 )
                        extended_score = score + step_q;
                    else
                        extended_score = score + kmer_size - (unsigned) std::floor(((double)dist * penalty));
                    if (extended_score > 0) {
                        if (!extending){ //first hit (2 kmers found at correct distance)
                            diag = curr_d.diag;
                            start = prev_d.wpos.pos;
                            end = curr_d.wpos.pos;
                            extending = true;
                            score = extended_score;
                        } else { //extending
                            end = curr_d.wpos.pos;
                            score = extended_score;
                        }
                    } else //stop extension score to join two kmers < 0
                        if (extending) {
                            if (end + kmer_size - start - 1 >= min_frag_size) {
                                count++;
                                frag.push_back(record_frag(start, end, diag,
                                                           score, numseqQ, curr_d.wpos.numSeq, count));
                            }
                            extending = false;
                            score=kmer_size;
                        }

                } else //stop extension diag or seq are different
                    if (extending) {
                        if (end + kmer_size - start - 1 >= min_frag_size) {
                            count++;
                            frag.push_back(record_frag(start, end, diag,
                                                       score, numseqQ, curr_d.wpos.numSeq, count));
                        }
                        extending = false;
                        score=kmer_size;
                }
                prev_d = curr_d;
            } //end size>2, diag_map loop
            if (extending) // Record hit at the end of the loop
                if (end + kmer_size - start - 1 >= min_frag_size) {
                    count++;
                    frag.push_back(record_frag(start, end, diag,
                                               score, numseqQ, curr_d.wpos.numSeq, count));
                }
        }

        if (verbose > 0) {
            std::cout << "Fragments number founds:" << count << std::endl;
        }
    }
}
//-------------------------------------------------------------------------
// Search for alignments with word matches
void Hasher::matchKmersHole(const BioSeq& sequence,
                            unsigned start, unsigned end, bool repeat,
                            std::vector< std::list<Diag> >& diag_map)
{
  unsigned last_pos=end-kmer_size;
  if(end<=kmer_size) return;

  std::string str=sequence.substr(start,end-start+1);
  const char* seq=str.c_str();

  unsigned key_d,dirhit=0;
  unsigned i=start;
  while(i<=last_pos) {
      bool found=false;
      key_d = hseq(seq);
      auto begin_d = hash2wpos[key_d];
      auto end_d = hash2wpos[key_d + 1];
      for (auto j = begin_d; j != end_d; j++) {
          if (j->numSeq == 0) continue;
          if (j->numSeq > 0) {
              long diag = long(i) - j->pos;
              if (!repeat || (repeat && i < j->pos)){
                  dirhit++;
                  diag_map[j->numSeq].push_back(Diag(diag, j->pos, j->numSeq));
                  found = true;
              }
          }
      }
      if(found){
      seq += step_q;
      i+=step_q;
      } else {
          seq += 1;
          i += 1;
      }
  }
  std::cout<<dirhit<<" hits found / ";
}
//-------------------------------------------------------------------------
// Search for alignments with word matches
void Hasher::matchKmersMinimizer(const BioSeq& sequence,
                            unsigned start, unsigned end, bool repeat,
                            std::vector< std::list<Diag> >& diag_map)
{
    unsigned last_pos=end-kmer_size;
    if(end<=kmer_size) return;

    std::string str=sequence.substr(start,end-start+1);
    const char* seq=str.c_str();

    unsigned key_d,dirhit=0;
    unsigned i=start;
    while(i<=last_pos) {
        bool found=false;
        key_d = mseq(seq);
        auto begin_d = hash2wpos[key_d];
        auto end_d = hash2wpos[key_d + 1];
        for (auto j = begin_d; j != end_d; j++) {
            if (j->numSeq == 0) continue;
            if (j->numSeq > 0) {
                long diag = long(i) - j->pos;
                if (!repeat || (repeat && i < j->pos)){
                    dirhit++;
                    diag_map[j->numSeq].push_back(Diag(diag, j->pos, j->numSeq));
                    found = true;
                }
            }
        }
        if(found){
            seq += step_q;
            i+=step_q;
        } else {
            seq += 1;
            i += 1;
        }
    }
    std::cout<<dirhit<<" hits found / ";
}
//-------------------------------------------------------------------------
// Search for Alignments
void Hasher::search(const BioSeq& sequence, unsigned start, unsigned end, unsigned numseq, unsigned connect_dist,
		unsigned min_frag_size, bool repeat, std::list< RangePair >& frag, unsigned verbose)
{
	clock_t clock_begin, clock_end;
	clock_begin = clock();
	std::cout<<"hashing query sequence..."<<std::flush;

	std::vector< std::list<Diag> > diag_map;
	diag_map.resize(subject_names.size()+1);
    if(algorithm==1)
        matchKmersHole(sequence, start, end, repeat, diag_map);
    else
        matchKmersMinimizer(sequence, start, end, repeat, diag_map);

	clock_end = clock();
	std::cout<<" --> Time spent: "<<(double)(clock_end-clock_begin)/CLOCKS_PER_SEC<<" seconds"<<std::endl;

	clock_begin = clock();
	std::cout<<"search fragments..."<<std::flush;
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
unsigned Hasher::fragCoverage(const std::list< RangePair >& frag)
{
    unsigned coverage=0;
    std::list< RangePair > frag_sort(frag);
    frag_sort.sort(RangePair::less);

    unsigned size=frag_sort.size();
    if(size>=2){
        auto curr_frag_it=frag_sort.begin();
        auto next_frag_it =curr_frag_it;
        next_frag_it++;
        while(next_frag_it != frag_sort.end()) {
            if (curr_frag_it->overlapQ(*next_frag_it)) {
                //TODO chose a subject name if different
                curr_frag_it->merge(*next_frag_it);
                next_frag_it = frag_sort.erase(next_frag_it);
            } else{
                curr_frag_it++;
                next_frag_it++;
            }
        }

    }
    for(auto & it : frag_sort)
        coverage+=it.getLength();
    return coverage;
}
//-------------------------------------------------------------------------
// Stats on rangePair lists
unsigned Hasher::fragScoreStat(const std::list< RangePair >& frag, double quantile, unsigned& coverage)
{
    coverage=0;
    if(frag.empty()){
        return 0;
    }
    std::vector<unsigned> score_list;
    for(const auto & curr_frag_it : frag) {
        unsigned len=curr_frag_it.getLength();
        coverage+=len;
        score_list.push_back(curr_frag_it.getScore());
    }
    sort(score_list.begin(), score_list.end());
    unsigned nb_frag=score_list.size();
    unsigned min_score=score_list.front();
    unsigned max_score=score_list.back();
    unsigned qval=score_list[(int)std::floor((double)score_list.size() * quantile)];
    std::cout << "Frag number=" << nb_frag << " / "
              << "min score=" << min_score << " / "
              << "max score=" << max_score << " / "
             <<"quantile ("<<quantile<<")="<<qval
             <<std::endl;
    return qval;
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
              <<"quantile ("<<quantile<<")="<<qval
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
void Hasher::fragScoreFilter(std::list< RangePair >& frag, unsigned min_score)
{
    std::cout<<"--Filter fragments score <"<<min_score<<" ... "<<std::flush;
    auto frag_it=frag.begin();
    while(frag_it != frag.end()) {
        if(frag_it->getScore()<min_score){
            frag_it = frag.erase(frag_it);
        }else{frag_it++;}
    }
    std::cout<<"done !"<<std::endl;
}
//-------------------------------------------------------------------------
// Set rangePair score and identity
void Hasher::fragSeqAlign(std::list< RangePair >& frag,
                          const SDGString& fasta_queryfilename, const SDGString& fasta_subjectfilename, bool reverse, unsigned verbose)
{
    FastaIstream query_in(fasta_queryfilename);
    if (!query_in) {
        std::cerr << "file:" << fasta_queryfilename << " does not exist!" << std::endl;
    }

    std::vector<BioSeq> subject_db;
    FastaIstream subject_in(fasta_subjectfilename);
    if (!query_in) {
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
                for(unsigned i=0;i<qlen;i++)
                {
                    if(qseq[i]==fragsseq[i]) count++;
                }

                curr_frag_it.setIdentity(((double)count)/qlen*100);
                curr_frag_it.setScore(count);
                if(verbose>0) std::cout << "Score = " << count<<" identity = " << ((double)count)/qlen * 100<< std::endl;
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
    char buff[2048];

    while (fileQ) {
        fileQ.getline(buff, 2048);
        if (*(buff) == '>') {
            SDGString chr_name(&buff[1]);
            chr_name = chr_name.trimL();
            chr_name = chr_name.trimR();
            num2nameQ.push_back(chr_name);
        }
    }
    fileQ.close();

    std::ifstream fileS(sfilename);
    while (fileS) {
        fileS.getline(buff, 2048);
        if (*(buff) == '>') {
            SDGString chr_name(&buff[1]);
            chr_name = chr_name.trimL();
            chr_name = chr_name.trimR();
            num2nameS.push_back(chr_name);
        }
    }
    fileS.close();

    for(auto & curr_frag_it : frag) {
        if(curr_frag_it.getNumQuery()>(long)num2nameQ.size()){
            std::cerr<<"Error query sequence number "<<curr_frag_it.getNumQuery()<<" doesn't exist!"<<std::endl;
            exit(EXIT_FAILURE);
        }
        SDGString qname=num2nameQ[curr_frag_it.getNumQuery()-1];

        if(curr_frag_it.getNumSubject()>(long)num2nameS.size()){
            std::cerr<<"Error subject sequence number "<<curr_frag_it.getNumSubject()<<" doesn't exist!"<<std::endl;
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
        for (const auto & curr_frag_it : frag) {
            if (curr_frag_it.getRangeQ().getNumChr() == numseq) {
                BioSeq sseq;
                if (curr_frag_it.getRangeQ().isPlusStrand()) {
                    sseq = seq.subseq(curr_frag_it.getRangeQ().getStart(), curr_frag_it.getRangeQ().getEnd() - curr_frag_it.getRangeQ().getStart() + 1);
                } else {
                    sseq = seq.subseq(curr_frag_it.getRangeQ().getEnd(), curr_frag_it.getRangeQ().getStart() - curr_frag_it.getRangeQ().getEnd() + 1);
                    sseq = sseq.complement();
                }
                std::istringstream subject_name(curr_frag_it.getRangeS().getNameSeq());
                std::string prefix_name;
                subject_name>>prefix_name;
                std::ostringstream name;
                name << prefix_name << " " <<seq.header << ":"
                    << curr_frag_it.getRangeQ().getStart() << ".."
                    << curr_frag_it.getRangeQ().getEnd();
                sseq.header=name.str();
                out << sseq;
            }
        }
    }
}
