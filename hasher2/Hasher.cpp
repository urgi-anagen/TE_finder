#include "Hasher.h"

//-------------------------------------------------------------------------
// Search for diagonal of word matches
void Hasher::diagSearch(const SDGBioSeq& sequence, std::vector< std::list<Diag> >& diag_map,std::vector< std::list<Diag> >& diag_map_comp,
  		unsigned connect_dist, unsigned kmer_size, unsigned min_frag_size, std::ostream& out, unsigned verbose)
{
  unsigned count=0;
  SDGString qname=sequence.getDE();

  unsigned len=sequence.length();
  for(std::vector< std::list<Diag> >::iterator iter_seq=diag_map.begin(); iter_seq!=diag_map.end(); iter_seq++)
  {
	  unsigned size=iter_seq->size();
	  if(size>2)
		{
		  iter_seq->sort();
		  SDGString sname=subject_db[iter_seq->front().wpos.numSeq-1].getDE();
		  unsigned start=0;
		  unsigned end=0;
		  unsigned score=0;
		  int diag=0;

		  std::list<Diag>::iterator iter_diag=iter_seq->begin();
		  Diag& prev_d=*iter_diag;
		  while( iter_diag!=iter_seq->end())
		  {
			  iter_diag++;
			  Diag& curr_d=*iter_diag;
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
						  score=1;
						}
					}
				  else //stop extension if distance between kmer too long
					  if(start!=0)
						{
						  if(end+kmer_size-start-1>=min_frag_size)
						  {
							 out<<qname<<"\t"
							 <<diag+start+1<<"\t"<<diag+end+kmer_size
							 <<"\t"
							 <<sname<<"\t"
							 <<start+1<<"\t"<<end+kmer_size
							 <<"\t0.0"
							 <<"\t"<<score
							 <<"\t"<<(double(score)*kmer_size/(end+kmer_size-start+1))*100
							 <<std::endl;
						  }
						  start=0;
						}
			  prev_d=curr_d;
		  } //end for
		  if(start!=0) // Record hit at the end of the loop
			{
			  if(end+kmer_size-start-1>=min_frag_size)
			  {
					 out<<qname<<"\t"
					 <<diag+start+1<<"\t"<<diag+end+kmer_size
					 <<"\t"
					 <<sname<<"\t"
					 <<start+1<<"\t"<<end+kmer_size
					 <<"\t0.0"
					 <<"\t"<<score
					 <<"\t"<<(double(score)*kmer_size/(end+kmer_size-start+1))*100
					 <<std::endl;
			  }
			}
		} //end size>2, diag_map loop
  }

  for(std::vector< std::list<Diag> >::iterator iter_seq=diag_map_comp.begin(); iter_seq!=diag_map_comp.end(); iter_seq++)
  {
	  unsigned size=iter_seq->size();
	  if(size>2)
		{
		  iter_seq->sort();
		  SDGString sname=subject_db[iter_seq->front().wpos.numSeq-1].getDE();

		  unsigned start=0;
		  unsigned end=0;
		  unsigned score=0;
		  int diag=0;

		  std::list<Diag>::iterator iter_diag=iter_seq->begin();
		  Diag& prev_d=*iter_diag;
		  while( iter_diag!=iter_seq->end())
		  {
			  iter_diag++;
			  Diag& curr_d=*iter_diag;
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
						  score=1;
						}
				  else //stop extension if distance between kmer too long
					  if(start!=0) // Record hit
						{
						  if(end+kmer_size-start-1>=min_frag_size)
						  {
							 out<<qname<<"\t"
							 <<len-(diag+start)<<"\t"<<len-(diag+end)-kmer_size+1
							 <<"\t"
							 <<sname<<"\t"
							 <<start+1<<"\t"<<end+kmer_size
							 <<"\t0.0"
							 <<"\t"<<score
							 <<"\t"<<(double(score)*kmer_size/(end+kmer_size-start+1))*100
							 <<std::endl;
						  }
						  start=0;
						}
			  prev_d=curr_d;
			} //end for
		  if(start!=0) // Record hit at the end of the loop
			{
			  if(end+kmer_size-start-1>=min_frag_size)
			  {
					 out<<qname<<"\t"
					 <<len-(diag+start)<<"\t"<<len-(diag+end)-kmer_size+1
					 <<"\t"
					 <<sname<<"\t"
					 <<start+1<<"\t"<<end+kmer_size
					 <<"\t0.0"
					 <<"\t"<<score
					 <<"\t"<<(double(score)*kmer_size/(end+kmer_size-start+1))*100
					 <<std::endl;
			  }
			}
		} //end size>2, diag_map_comp loop

  }
  if(verbose>0)
  {
	  std::cout<<"Fragments number founds:"<<count<<std::endl;
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

  unsigned key_d=0,key_c=0;
  for(unsigned i=start;i<=last_pos;i+=step_q)
    {
      key_d=hseq(seq);
      seq+=step_q;

      std::vector<KmerSpos>::iterator begin_d=hash2wpos[key_d];
      std::vector<KmerSpos>::iterator end_d=hash2wpos[key_d+1];
      for(std::vector<KmerSpos>::iterator j=begin_d;j!=end_d;j++)
	  diag_map[j->numSeq].push_back(Diag(i-j->pos,j->pos,j->numSeq));

      key_c=hseq(seq_comp);
      seq_comp+=step_q;

      std::vector<KmerSpos>::iterator begin_c=hash2wpos[key_c];
      std::vector<KmerSpos>::iterator end_c=hash2wpos[key_c+1];
      for(std::vector<KmerSpos>::iterator j=begin_c;j!=end_c;j++)
	  diag_map_comp[j->numSeq].push_back(Diag(i-j->pos,j->pos,j->numSeq));

    }
}
//-------------------------------------------------------------------------
void Hasher::search(const SDGBioSeq& sequence, unsigned start, unsigned end, unsigned numseq, unsigned connect_dist,
		unsigned min_frag_size, bool repeat, std::ostream& out, unsigned verbose)
{
	clock_t clock_begin, clock_end;
	clock_begin = clock();
	std::cout<<"hashing query sequence..."<<std::flush;

	std::vector< std::list<Diag> > diag_map, diag_map_comp;
	diag_map.resize(subject_db.size()+1);
	diag_map_comp.resize(subject_db.size()+1);
	matchKmers(sequence, start, end, numseq, repeat, diag_map, diag_map_comp);
	std::cout<<"ok"<<std::endl;
	std::cout<<diag_map.size()<<" direct hits found / ";
	std::cout<<diag_map_comp.size()<<" reverse hits found";

	clock_end = clock();
	std::cout<<" --> Time spent: "<<(double)(clock_end-clock_begin)/CLOCKS_PER_SEC<<" seconds"<<std::endl;

	clock_begin = clock();
	std::cout<<"search fragments..."<<std::flush;
	diagSearch(sequence, diag_map, diag_map_comp, connect_dist, kmer_size, min_frag_size, out, verbose);
	diag_map.clear();
	std::cout<<"ok"<<std::endl;
	clock_end = clock();
	std::cout<<" --> Time spent: "<<(double)(clock_end-clock_begin)/CLOCKS_PER_SEC<<" seconds"<<std::endl;
}




