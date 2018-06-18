#include "Hasher.h"

//-------------------------------------------------------------------------
// Search for diagonal of word matches
void Hasher::diagSearch(const SDGBioSeq& sequence, std::vector< Diag >& diag_map,std::vector< Diag >& diag_map_comp,
  		unsigned connect_dist, unsigned kmer_size, unsigned min_frag_size, unsigned verbose)
{
  frag.clear();
  unsigned count=0;
  unsigned size=diag_map.size();
  unsigned len=sequence.length();
  if(size>2)
    {
      sort(diag_map.begin(),diag_map.end());

      unsigned start=0;
      unsigned end=0;
      unsigned numSeq=0;
      unsigned score=0;
      int diag=0;

      Diag& prev_d=diag_map[0];
      for( unsigned i=1; i<size; ++i)
      {
		  Diag& curr_d=diag_map[i];
		  if(prev_d.diag==curr_d.diag
			 && prev_d.wpos.numSeq==curr_d.wpos.numSeq
			 && prev_d.wpos.pos+connect_dist>=curr_d.wpos.pos)
			  	{
				  if(start!=0) //extending
					{
					  end=curr_d.wpos.pos;
					  score++;
					}
				  else //first hit (2 kmers found at correct distance)
					{
					  diag=prev_d.diag;
					  start=prev_d.wpos.pos;
					  end=curr_d.wpos.pos;
					  numSeq=prev_d.wpos.numSeq;
					  score=1;
					}
			  	}
			  else //stop extension if distance between kmer too long
				  if(start!=0)
					{
					  if(end+kmer_size-start-1>=min_frag_size)
					  {
						  RangeAlign r1(nbseqQ,diag+start+1,diag+end+kmer_size);
						  RangeAlign r2(numSeq,start+1,end+kmer_size);

						  RangePair rp(r1,r2);
						  rp.setId(++count);
						  rp.setScore(score);
						  rp.setIdentity(double(score)*kmer_size/r1.getLength());
						  rp.setLength(r1.getLength());
						  frag.push_back(rp);
					  }
					  start=0;
					}
		  prev_d=curr_d;
      } //end for
      if(start!=0) // Record hit at the end of the loop
		{
		  if(end+kmer_size-start-1>=min_frag_size)
		  {
			  RangeAlign r1(nbseqQ,diag+start+1,diag+end+kmer_size);
			  RangeAlign r2(numSeq,start+1,end+kmer_size);

			  RangePair rp(r1,r2);
			  rp.setId(++count);
			  rp.setScore(score);
			  rp.setIdentity(double(score)*kmer_size/r1.getLength());
			  rp.setLength(r1.getLength());
			  frag.push_back(rp);
		  }
		}
    } //end size>2, diag_map loop

  size=diag_map_comp.size();
  if(size>2)
    {
      sort(diag_map_comp.begin(),diag_map_comp.end());

      unsigned start=0;
      unsigned end=0;
      unsigned numSeq=0;
      unsigned score=0;
      int diag=0;

      Diag& prev_d=diag_map_comp[0];
      for( unsigned i=1; i<size; ++i)
		{
		  Diag& curr_d=diag_map_comp[i];
		  if(prev_d.diag==curr_d.diag
			 && prev_d.wpos.numSeq==curr_d.wpos.numSeq
			  && prev_d.wpos.pos+connect_dist>=curr_d.wpos.pos)
				  if(start!=0) //extending
					{
					  end=curr_d.wpos.pos;
					  score++;
					}
				  else //first hit (2 kmers found at correct distance)
					{
					  diag=prev_d.diag;
					  start=prev_d.wpos.pos;
					  end=curr_d.wpos.pos;
					  numSeq=prev_d.wpos.numSeq;
					  score=1;
					}
			  else //stop extension if distance between kmer too long
				  if(start!=0) // Record hit
					{
					  if(end+kmer_size-start-1>=min_frag_size)
					  {
						  RangeAlign r1(nbseqQ,len-(diag+start),
								len-(diag+end)-kmer_size+1);
						  RangeAlign r2(numSeq,start+1,end+kmer_size);

						  RangePair rp(r1,r2);
						  rp.setId(++count);
						  rp.setScore(score);
						  rp.setIdentity(double(score)*kmer_size/r1.getLength());
						  rp.setLength(r1.getLength());

						  frag.push_back(rp);
					  }
					  start=0;
					}
		  prev_d=curr_d;
		} //end for
      if(start!=0) // Record hit at the end of the loop
		{
		  if(end+kmer_size-start-1>=min_frag_size)
		  {
			  RangeAlign r1(nbseqQ,len-(diag+start),
						len-(diag+end)-kmer_size+1);
			  RangeAlign r2(numSeq,start+1,end+kmer_size);


			  RangePair rp(r1,r2);
			  rp.setId(++count);
			  rp.setScore(score);
			  rp.setIdentity(double(score)*kmer_size/r1.getLength());
			  rp.setLength(r1.getLength());

			  frag.push_back(rp);
		  }
		}
    } //end size>2, diag_map_comp loop


  if(verbose>0)
  {
	  std::cout<<"Fragments number founds:"<<count<<std::endl;
  }
}
//-------------------------------------------------------------------------
// Search for alignments with word matches
void Hasher::matchKmers(const SDGBioSeq& sequence,
	    unsigned start, unsigned end, unsigned numseq, bool repeat,
		std::vector< Diag >& diag_map, std::vector< Diag >& diag_map_comp)
{
  unsigned last_pos=end-kmer_size;
  if(end<=kmer_size) return;

  std::string str=sequence.toString().substr(start,end-start);
  const char* seq=str.c_str();

  SDGBioSeq comp_sequence=newSDGMemBioSeq(str).complement();
  std::string str_comp=comp_sequence.toString();
  const char* seq_comp=str_comp.c_str();

  unsigned key_d=0,key_c=0;
  for(unsigned i=start;i<=last_pos;i+=step_q)
    {
      key_d=hseq(seq);
      seq+=step_q;

      std::vector<KmerSpos>::iterator begin_d=hash2wpos[key_d];
      std::vector<KmerSpos>::iterator end_d=hash2wpos[key_d+1];
      for(std::vector<KmerSpos>::iterator j=begin_d;j!=end_d;j++)
	  diag_map.push_back(Diag(i-j->pos,j->pos,j->numSeq));

      key_c=hseq(seq_comp);
      seq_comp+=step_q;

      std::vector<KmerSpos>::iterator begin_c=hash2wpos[key_c];
      std::vector<KmerSpos>::iterator end_c=hash2wpos[key_c+1];
      for(std::vector<KmerSpos>::iterator j=begin_c;j!=end_c;j++)
	  diag_map_comp.push_back(Diag(i-j->pos,j->pos,j->numSeq));

    }
}
//-------------------------------------------------------------------------
void Hasher::search(const SDGBioSeq& sequence, unsigned start, unsigned end, unsigned numseq, unsigned connect_dist,
		unsigned min_frag_size, bool repeat, unsigned verbose)
{
	clock_t clock_begin, clock_end;
	clock_begin = clock();
	std::cout<<"hashing query sequence..."<<std::flush;

	std::vector< Diag > diag_map, diag_map_comp;
	matchKmers(sequence, start, end, numseq, repeat, diag_map, diag_map_comp);
	std::cout<<"ok"<<std::endl;
	std::cout<<diag_map.size()<<" direct hits found / ";
	std::cout<<diag_map_comp.size()<<" reverse hits found";

	clock_end = clock();
	std::cout<<" --> Time spent: "<<(double)(clock_end-clock_begin)/CLOCKS_PER_SEC<<" seconds"<<std::endl;

	clock_begin = clock();
	std::cout<<"search fragments..."<<std::flush;
	diagSearch(sequence, diag_map, diag_map_comp, connect_dist, kmer_size, min_frag_size, verbose);
	diag_map.clear();
	std::cout<<"ok"<<std::endl;
	std::cout<<frag.size()<<" fragments found";
	clock_end = clock();
	std::cout<<" --> Time spent: "<<(double)(clock_end-clock_begin)/CLOCKS_PER_SEC<<" seconds"<<std::endl;
}
//-------------------------------------------------------------------------
// Write aligned fragments in .align format
void Hasher::write_align(const SDGBioSeq& sequence, std::ostream& out)
{
  SDGString qname=sequence.getDE();
  for(std::list< RangePair >::iterator i=frag.begin();i!=frag.end();
      i++)
    {
		  SDGString sname=subject_db[i->getRangeS().getNumChr()-1].getDE();
		 out<<qname<<"\t"
		 <<i->getRangeQ().getStart()<<"\t"<<i->getRangeQ().getEnd()
		 <<"\t"
		 <<sname<<"\t"
		 <<i->getRangeS().getStart()<<"\t"<<i->getRangeS().getEnd()
		 <<"\t"<<i->getE_value()
		 <<"\t"<<i->getScore()
		 <<"\t"<<i->getIdentity()
		 <<std::endl;
    }
}





