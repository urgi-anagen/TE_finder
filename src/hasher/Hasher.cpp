#include "Hasher.h"

//-------------------------------------------------------------------------
// Search for diagonal of word matches
void Hasher::diagSearch(const SDGBioSeq &sequence, unsigned numseqQ, std::vector<std::list<Diag> > &diag_map,
                        std::vector<std::list<Diag> > &diag_map_comp,
                        unsigned connect_dist, unsigned kmer_size, unsigned min_frag_size,
                        std::list< RangePair >& frag, unsigned verbose) {
    unsigned count = 0;
    SDGString qname = sequence.getDE();

    unsigned len = sequence.length();
    unsigned curr_seq=0;
    for (std::vector<std::list<Diag> >::iterator iter_seq = diag_map.begin(); iter_seq != diag_map.end(); iter_seq++) {
        unsigned size = iter_seq->size();
        if (size > 2) {
            iter_seq->sort();
            unsigned start = 0;
            unsigned end = 0;
            unsigned score = 0;
            int diag = 0;

            std::list<Diag>::iterator iter_diag = iter_seq->begin();
            Diag prev_d = *iter_diag;
            while (++iter_diag != iter_seq->end()) {
                Diag curr_d = *iter_diag;
                curr_seq=curr_d.wpos.numSeq;

                if (prev_d.diag == curr_d.diag
                    && prev_d.wpos.numSeq == curr_d.wpos.numSeq
                    && prev_d.wpos.pos + connect_dist >= curr_d.wpos.pos) {
                    if (start != 0) //extending
                    {
                        end = curr_d.wpos.pos;
                        score++;
                    } else //first hit (2 kmers found at correct distance)
                    {
                        diag = prev_d.diag;
                        start = prev_d.wpos.pos;
                        end = curr_d.wpos.pos;
                        score = 1;
                    }
                } else //stop extension if distance between kmer too long
                if (start != 0) {
                    if (end + kmer_size - start - 1 >= min_frag_size) {
                        count++;

                        unsigned qstart = diag + start + 1;
                        unsigned qend = diag + end + kmer_size;
                        unsigned sstart = start + 1;
                        unsigned send = end + kmer_size;

                        RangePair rp= rangePairFactory(numseqQ, qstart, qend, curr_seq , sstart, send, score, step_q, count);
                        frag.push_back(rp);

                    }
                    start = 0;
                }
                prev_d = curr_d;
            } //end for
            if (start != 0) // Record hit at the end of the loop
            {
                if (end + kmer_size - start - 1 >= min_frag_size) {
                    count++;
                    unsigned qstart = diag + start + 1;
                    unsigned qend = diag + end + kmer_size;
                    unsigned sstart = start + 1;
                    unsigned send = end + kmer_size;
                    RangePair rp= rangePairFactory(numseqQ, qstart, qend, curr_seq , sstart, send, score, step_q, count);
                    frag.push_back(rp);
                }
            }
        } //end size>2, diag_map loop
    }

    curr_seq=0;
    for (std::vector<std::list<Diag> >::iterator iter_seq = diag_map_comp.begin();
         iter_seq != diag_map_comp.end(); iter_seq++) {
        unsigned size = iter_seq->size();
        if (size > 2) {
            iter_seq->sort();

            unsigned start = 0;
            unsigned end = 0;
            unsigned score = 0;
            int diag = 0;

            std::list<Diag>::iterator iter_diag = iter_seq->begin();
            Diag prev_d = *iter_diag;
            while (++iter_diag != iter_seq->end()) {
                Diag curr_d = *iter_diag;
                curr_seq=curr_d.wpos.numSeq;

                if (prev_d.diag == curr_d.diag
                    && prev_d.wpos.numSeq == curr_d.wpos.numSeq
                    && prev_d.wpos.pos + connect_dist >= curr_d.wpos.pos)
                    if (start != 0) //extending
                    {
                        end = curr_d.wpos.pos;
                        score++;
                    } else //first hit (2 kmers found at correct distance)
                    {
                        diag = prev_d.diag;
                        start = prev_d.wpos.pos;
                        end = curr_d.wpos.pos;
                        score = 1;
                    }
                else //stop extension if distance between kmer too long
                if (start != 0) // Record hit
                {
                    if (end + kmer_size - start - 1 >= min_frag_size) {
                        count++;
                        unsigned qstart = len - (diag + start);
                        unsigned qend = len - (diag + end) - kmer_size + 1;
                        unsigned sstart = start + 1;
                        unsigned send = end + kmer_size;
                        RangePair rp= rangePairFactory(numseqQ, qstart, qend, curr_seq , sstart, send, score, step_q, count);
                        frag.push_back(rp);
                    }
                    start = 0;
                }
                prev_d = curr_d;
            } //end for
            if (start != 0) // Record hit at the end of the loop
            {
                if (end + kmer_size - start - 1 >= min_frag_size) {
                    count++;
                    unsigned qstart = len - (diag + start);
                    unsigned qend = len - (diag + end) - kmer_size + 1;
                    unsigned sstart = start + 1;
                    unsigned send = end + kmer_size;
                    RangePair rp= rangePairFactory(numseqQ, qstart, qend, curr_seq , sstart, send, score, step_q, count);
                    frag.push_back(rp);
                }
            }
        } //end size>2, diag_map_comp loop

    }
    if (verbose > 0) {
        std::cout << "Fragments number founds:" << count << std::endl;
    }
}




//-------------------------------------------------------------------------
// Search for alignments with word matches
void Hasher::matchKmers(const SDGBioSeq& sequence,
	    unsigned start, unsigned end, unsigned numseq, bool repeat,
		std::vector< std::list<Diag> >& diag_map, std::vector< std::list<Diag> >& diag_map_comp)
{
  unsigned last_pos=end-kmer_size;
  if(end<=kmer_size) return;

  std::string str=sequence.toString().substr(start,end-start);
  const char* seq=str.c_str();

  SDGBioSeq comp_sequence=newSDGMemBioSeq(str).complement();
  std::string str_comp=comp_sequence.toString();
  const char* seq_comp=str_comp.c_str();

  unsigned key_d=0,dirhit=0;
  unsigned i=start;
  while(i<=last_pos) {
      bool found=false;
      key_d = hseq(seq);
      std::vector<KmerSpos>::iterator begin_d = hash2wpos[key_d];
      std::vector<KmerSpos>::iterator end_d = hash2wpos[key_d + 1];
      for (std::vector<KmerSpos>::iterator j = begin_d; j != end_d; j++) {
          if (j->numSeq > 0) {
              dirhit++;
              long diag = long(i) - j->pos;
              diag_map[j->numSeq].push_back(Diag(diag, j->pos, j->numSeq));
              found=true;
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
  std::cout<<dirhit<<" direct hits found / ";

  unsigned key_c=0,comphit=0;
  i=start;
  while(i<=last_pos) {
      bool found=false;
      key_c=hseq(seq_comp);
      std::vector<KmerSpos>::iterator begin_c=hash2wpos[key_c];
      std::vector<KmerSpos>::iterator end_c=hash2wpos[key_c+1];
      for(std::vector<KmerSpos>::iterator j=begin_c;j!=end_c;j++)
      {
          if (j->numSeq > 0) {
              comphit++;
              long diag = long(i) - j->pos;
              diag_map_comp[j->numSeq].push_back(Diag(diag, j->pos, j->numSeq));
              found=true;
          }
      }
      if(found){
          seq_comp += step_q;
          i+=step_q;
      } else {
          seq_comp += 1;
          i += 1;
      }
    }

	std::cout<<comphit<<" reverse hits found"<<std::endl;;
}
//-------------------------------------------------------------------------
// Search for Alignments
void Hasher::search(const SDGBioSeq& sequence, unsigned start, unsigned end, unsigned numseq, unsigned connect_dist,
		unsigned min_frag_size, bool repeat, std::list< RangePair >& frag, unsigned verbose)
{
	clock_t clock_begin, clock_end;
	clock_begin = clock();
	std::cout<<"hashing query sequence..."<<std::flush;

	std::vector< std::list<Diag> > diag_map, diag_map_comp;
	diag_map.resize(subject_names.size()+1);
	diag_map_comp.resize(subject_names.size()+1);
	matchKmers(sequence, start, end, numseq, repeat, diag_map, diag_map_comp);

	clock_end = clock();
	std::cout<<" --> Time spent: "<<(double)(clock_end-clock_begin)/CLOCKS_PER_SEC<<" seconds"<<std::endl;

	clock_begin = clock();
	std::cout<<"search fragments..."<<std::flush;
	diagSearch(sequence, numseq, diag_map, diag_map_comp, connect_dist, kmer_size, min_frag_size, frag, verbose);
	diag_map.clear();
	std::cout<<"ok"<<std::endl;
	clock_end = clock();
	std::cout<<" --> Time spent: "<<(double)(clock_end-clock_begin)/CLOCKS_PER_SEC<<" seconds"<<std::endl;

    clock_begin = clock();
    std::cout<<"merge fragments..."<<std::flush;
    fragMerge(frag);
    std::cout<<"ok"<<std::endl;
    std::cout<<frag.size()<<" ranges found";
    clock_end = clock();
    std::cout<<" --> Time spent: "<<(double)(clock_end-clock_begin)/CLOCKS_PER_SEC<<" seconds"<<std::endl;
}
//-------------------------------------------------------------------------
// Merge overlapping rangePair
void Hasher::fragMerge(std::list< RangePair >& frag)
{

    frag.sort(RangePair::less);
    unsigned size=frag.size();

    if(size>=2){
        std::list< RangePair >::iterator curr_frag_it=frag.begin();
        std::list< RangePair >::iterator next_frag_it =curr_frag_it;
        next_frag_it++;
        while( next_frag_it != frag.end()) {

            if (curr_frag_it->overlapQ(*next_frag_it)) {
                curr_frag_it->merge(*next_frag_it);
                next_frag_it = frag.erase(next_frag_it);
            } else{
                curr_frag_it++;
                next_frag_it++;
            }
        }

    }
}
//-------------------------------------------------------------------------
// Stats on rangePair lists
unsigned Hasher::fragStat(const std::list< RangePair >& frag, double quantile, unsigned& coverage)
{

    coverage=0;
    if(frag.size()==0){
        return 0;
    }
    std::vector<unsigned> length_list;
    for( std::list< RangePair >::const_iterator curr_frag_it=frag.begin();
         curr_frag_it != frag.end() ; curr_frag_it++) {
        unsigned len=curr_frag_it->getLength();
        coverage+=len;
        length_list.push_back(len);
    }
    sort(length_list.begin(),length_list.end());
    unsigned nb_frag=length_list.size();
    unsigned min_len=length_list.front();
    unsigned max_len=length_list.back();
    unsigned qval=length_list[(int)std::floor(length_list.size()*quantile)];
    std::cout<<"Frag number="<<nb_frag<<" / "
             <<"min length="<<min_len<<" / "
             <<"max length="<<max_len<<" / "
             <<"quantile ("<<quantile<<")="<<qval
             <<std::endl;
    return qval;
}
//-------------------------------------------------------------------------
// Stats on rangePair lists
void Hasher::fragLenFilter(std::list< RangePair >& frag, unsigned min_len)
{
    std::cout<<"--Filter fragments <"<<min_len<<" ... "<<std::flush;
    std::list< RangePair >::iterator frag_it=frag.begin();
    while(frag_it != frag.end()) {
        if(frag_it->getLength()<min_len){
            frag_it = frag.erase(frag_it);
        }else{frag_it++;}
    }
    std::cout<<"done !"<<std::endl;
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

    for( std::list< RangePair >::iterator curr_frag_it=frag.begin();
         curr_frag_it != frag.end() ; curr_frag_it++) {
        if(curr_frag_it->getNumQuery()>num2nameQ.size()){
            std::cerr<<"Error query sequence number "<<curr_frag_it->getNumQuery()<<" doesn't exist!"<<std::endl;
            exit(EXIT_FAILURE);
        }
        SDGString qname=num2nameQ[curr_frag_it->getNumQuery()-1];

        if(curr_frag_it->getNumSubject()>num2nameS.size()){
            std::cerr<<"Error subject sequence number "<<curr_frag_it->getNumSubject()<<" doesn't exist!"<<std::endl;
            exit(EXIT_FAILURE);
        }
        SDGString sname=num2nameS[curr_frag_it->getNumSubject()-1];;
        curr_frag_it->setQSName(qname,sname);
        curr_frag_it->writetxt(out);
    }
}
//-------------------------------------------------------------------------
// Write a rangePair lists
void Hasher::fragSeqWrite(const std::list< RangePair >& frag, const SDGString& fasta_filename, SDGFastaOstream& out)
{
    SDGFastaIstream in(fasta_filename);
    if (!in) {
        std::cerr << "file:" << fasta_filename << " does not exist!" << std::endl;
    }
    unsigned numseq=0;
    while (in) {
        SDGBioSeq seq;
        if (in)
            in >> seq;
        numseq++;
        std::cout << seq.getDE() << " len:" << seq.length() << " read!" << std::endl;
        for (std::list<RangePair>::const_iterator curr_frag_it = frag.begin();
             curr_frag_it != frag.end(); curr_frag_it++) {
            if (curr_frag_it->getRangeQ().getNumChr() == numseq) {
                SDGBioSeq sseq;
                if (curr_frag_it->getRangeQ().isPlusStrand()) {
                    sseq = seq.subseq(curr_frag_it->getRangeQ().getStart(), curr_frag_it->getRangeQ().getEnd() - curr_frag_it->getRangeQ().getStart() + 1);
                } else {
                    sseq = seq.subseq(curr_frag_it->getRangeQ().getEnd(), curr_frag_it->getRangeQ().getStart() - curr_frag_it->getRangeQ().getEnd() + 1);
                    sseq = sseq.complement();
                }
                std::istringstream subject_name(curr_frag_it->getRangeS().getNameSeq());
                std::string prefix_name;
                subject_name>>prefix_name;
                std::ostringstream name;
                name << prefix_name << " " <<seq.getDE() << ":"
                    << curr_frag_it->getRangeQ().getStart() << ".."
                    << curr_frag_it->getRangeQ().getEnd();
                sseq.setDE((SDGString) name.str());
                out << sseq;
            }
        }
    }
}
