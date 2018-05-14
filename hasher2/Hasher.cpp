#include "Hasher.h"

//-------------------------------------------------------------------------
// Search for diagonal of word matches
void Hasher::diagSearch(const SDGBioSeq& sequence, std::vector< Diag >& diag_map,std::vector< Diag >& diag_map_comp,
  		unsigned connect_dist, unsigned kmer_size)
{
	std::vector< std::pair<unsigned,RangePair> > frag;
  unsigned count=0;
  unsigned size=diag_map.size();
  unsigned len=sequence.length();
  unsigned wdlen=wdist*kmer_size;
  if(size>2)
    {
      sort(diag_map.begin(),diag_map.end());

      unsigned start=0;
      unsigned end=0;
      unsigned numSeq=0;
      int diag=0;
      bool in=false;
      bool found=false;

      Diag& prev_d=diag_map[0];
      for( unsigned i=1; i<size; ++i)
	{
	  Diag& curr_d=diag_map[i];
	  if(prev_d.diag==curr_d.diag
	     && prev_d.wpos.numSeq==curr_d.wpos.numSeq)
	    {
	      if(prev_d.wpos.pos+wdlen>=curr_d.wpos.pos)
		{
		  if(in)
		    {
		      end=curr_d.wpos.pos;
		    }
		  else
		    {
		      diag=prev_d.diag;
		      start=prev_d.wpos.pos;
		      end=curr_d.wpos.pos;
		      numSeq=prev_d.wpos.numSeq;
		      in=true;
		    }
		}
	      else
		{
		  if(in) found=true;
		  in=false;
		}
	    }
	  else
	    {
	      if(in) found=true;
	      in=false;
	    }
	  if(found)
	    {
	      RangeAlign r1(nbseqQ,diag+start+1,diag+end+kmer_size);
	      RangeAlign r2(numSeq,start+1,end+kmer_size);

	      RangePair rp(r1,r2);
	      rp.setId(++count);
	      rp.setScore(r1.getLength());
	      rp.setIdentity(1.00);
	      rp.setLength(r1.getLength());
	      frag.push_back(std::pair<unsigned,RangePair>(rp.getScore(),rp));
	      found=false;
	    }
	  prev_d=curr_d;
	}
      if(in)
	{
	  RangeAlign r1(nbseqQ,diag+start+1,diag+end+kmer_size);
	  RangeAlign r2(numSeq,start+1,end+kmer_size);

	  RangePair rp(r1,r2);
	  rp.setId(++count);
	  rp.setScore(r1.getLength());
	  rp.setIdentity(1.00);
	  rp.setLength(r1.getLength());
	  frag.push_back(std::pair<unsigned,RangePair>(rp.getScore(),rp));
	}
    }

  size=diag_map_comp.size();
  if(size>2)
    {
      sort(diag_map_comp.begin(),diag_map_comp.end());

      unsigned start=0;
      unsigned end=0;
      unsigned numSeq=0;
      int diag=0;
      bool in=false;
      bool found=false;
      Diag& prev_d=diag_map_comp[0];
      for( unsigned i=1; i<size; ++i)
	{
	  Diag& curr_d=diag_map_comp[i];
	  if(prev_d.diag==curr_d.diag
	     && prev_d.wpos.numSeq==curr_d.wpos.numSeq)
	    {
	      if(prev_d.wpos.pos+wdlen>=curr_d.wpos.pos)
		{
		  if(in)
		    {
		      end=curr_d.wpos.pos;
		    }
		  else
		    {
		      diag=prev_d.diag;
		      start=prev_d.wpos.pos;
		      end=curr_d.wpos.pos;
		      numSeq=prev_d.wpos.numSeq;
		      in=true;
		    }
		}
	      else
		{
		  if(in) found=true;
		  in=false;
		}
	    }
	  else
	    {
	      if(in) found=true;
	      in=false;
	    }
	  if(found)
	    {
	      RangeAlign r1(nbseqQ,len-(diag+start),
			    len-(diag+end)-kmer_size+1);
	      RangeAlign r2(numSeq,start+1,end+kmer_size);

	      RangePair rp(r1,r2);
	      rp.setId(++count);
	      rp.setScore(r1.getLength());
	      rp.setIdentity(1.00);
	      rp.setLength(r1.getLength());

	      frag.push_back(std::pair<unsigned,RangePair>(rp.getScore(),rp));
	      found=false;
	    }
	  prev_d=curr_d;
	}
      if(in)
	{
	  RangeAlign r1(nbseqQ,len-(diag+start),
			    len-(diag+end)-kmer_size+1);
	  RangeAlign r2(numSeq,start+1,end+kmer_size);


	  RangePair rp(r1,r2);
	  rp.setId(++count);
	  rp.setScore(r1.getLength());
	  rp.setIdentity(1.00);
	  rp.setLength(r1.getLength());

	  frag.push_back(std::pair<unsigned,RangePair>(rp.getScore(),rp));
	}
    }

  unsigned count2=0;
  map_align.clear();

  sort(frag.begin(),frag.end(),comp());

  for(std::vector< std::pair<unsigned,RangePair> >::reverse_iterator i=frag.rbegin();
      i!=frag.rend() && count2<nfrag;i++)
    {
      count2++;
      RangePair rangePair(i->second);
      std::list<RangePair>& al_list
	=map_align[Key(rangePair.getRangeQ().getNumChr(),
		       rangePair.getRangeS().getNumChr())];

      al_list.insert(al_list.begin(),rangePair);
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
void Hasher::search(const SDGBioSeq& sequence, unsigned start, unsigned end, unsigned numseq, bool repeat)
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
	diagSearch(sequence, diag_map, diag_map_comp,(wdist+1)*kmer_size,kmer_size);
	diag_map.clear();
	std::cout<<"ok"<<std::endl;
	std::cout<<map_align.size()<<" fragments found";
	clock_end = clock();
	std::cout<<" --> Time spent: "<<(double)(clock_end-clock_begin)/CLOCKS_PER_SEC<<" seconds"<<std::endl;
}
//-------------------------------------------------------------------------
void Hasher::fragAlign(double match,double mism, double gopen,
			   double gext, unsigned over, bool join)
{
  std::cout<<"fragment connexion ..."<<std::flush;
  map_path.clear();
  if(join)
    {
      for(MapAlign::iterator m=map_align.begin(); m!=map_align.end();m++)
		{
		  FragAlign fragAlign(mism,gopen,gext,over);
		  map_path[m->first]=fragAlign.join(m->second);
		  m->second.clear();
		}
    }
  else
    {
      for(MapAlign::iterator m=map_align.begin(); m!=map_align.end();m++)
		{
		  std::list<RangePairSet> path;
		  for(std::list<RangePair>::iterator i=m->second.begin();
			  i!=m->second.end();i++)
			{
			  path.push_back(RangePairSet(*i));
			}
		  map_path[m->first]=path;
		  m->second.clear();
		}
    }
  map_align.clear();
  std::cout<<"ok"<<std::endl;
}
//-------------------------------------------------------------------------
void Hasher::print_frag(const SDGBioSeq& sequence, std::ostream& out)
{
  SDGString qname=sequence.getDE();
  for(MapAlign::iterator iter_hash=map_align.begin();iter_hash!=map_align.end();
      iter_hash++)
    {
      SDGString sname=subject_db[iter_hash->first.second-1].getDE();
      for(std::list<RangePair>::iterator i=iter_hash->second.begin();
	  i!=iter_hash->second.end();i++)
		{
		  out<<qname<<"\t"
			 <<i->getRangeQ().getStart()<<".."<<i->getRangeQ().getEnd()
			 <<"\t"
			 <<sname<<"\t"
			 <<i->getRangeS().getStart()<<".."<<i->getRangeS().getEnd()
			 <<"\t"
			 <<i->getScore()<<std::endl;
		}
    }
}
//-------------------------------------------------------------------------
void Hasher::print(const SDGBioSeq& sequence, unsigned min_size,std::ostream& out)
{
  SDGString qname=sequence.getDE();

  for(MapPath::iterator iter_hash=map_path.begin();iter_hash!=map_path.end();
      iter_hash++)
    {
      SDGString sname=subject_db[iter_hash->first.second-1].getDE();
      for(std::list<RangePairSet>::iterator i=iter_hash->second.begin();
	  i!=iter_hash->second.end();i++)
		{
		  if(i->getRangeQ().getLength()>min_size
			 && i->getRangeS().getLength()>min_size)
			{
			  out<<qname<<"\t"
			 <<i->getRangeQ().getStart()<<".."<<i->getRangeQ().getEnd()
			 <<"\t"
			 <<sname<<"\t"
			 <<i->getRangeS().getStart()<<".."<<i->getRangeS().getEnd()
			 <<"\t"
			 <<i->getScore()<<std::endl;
			}
		}
    }
}
//-------------------------------------------------------------------------
void Hasher::write(const SDGBioSeq& sequence, unsigned min_size,std::ostream& out)
{
  SDGString qname=sequence.getDE();
  for(MapPath::iterator iter_hash=map_path.begin();iter_hash!=map_path.end();
      iter_hash++)
    {
      SDGString sname=subject_db[iter_hash->first.second-1].getDE();
      for(std::list<RangePairSet>::iterator i=iter_hash->second.begin();
	  i!=iter_hash->second.end();i++)
		{
		  if(i->getRangeQ().getLength()>min_size
			 && i->getRangeS().getLength()>min_size)
			{
			  out<<qname<<"\t"
			 <<i->getRangeQ().getStart()<<"\t"<<i->getRangeQ().getEnd()
			 <<"\t"
			 <<sname<<"\t"
			 <<i->getRangeS().getStart()<<"\t"<<i->getRangeS().getEnd()
			 <<"\t"
			 <<i->getScore()<<std::endl;
			}
		}
    }
}
//-------------------------------------------------------------------------
void Hasher::write_align(const SDGBioSeq& sequence, unsigned min_size,std::ostream& out)
{
  SDGString qname=sequence.getDE();
  for(MapPath::iterator iter_hash=map_path.begin();iter_hash!=map_path.end();
      iter_hash++)
    {
      SDGString sname=subject_db[iter_hash->first.second-1].getDE();
      for(std::list<RangePairSet>::iterator i=iter_hash->second.begin();
	  i!=iter_hash->second.end();i++)
		{
		  if(i->getRangeQ().getLength()>min_size
			 && i->getRangeS().getLength()>min_size)
			{
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
    }
}
//-------------------------------------------------------------------------
void Hasher::extend(const SDGBioSeq& sequence, const SDGBioSeq& comp_sequence, unsigned min_size, unsigned verbose)
{
  FastExtAlign fastExtAlign;
  FastExtAlign rfastExtAlign;
  unsigned seqlen=sequence.length();

  fastExtAlign.setSeq1(sequence);
  rfastExtAlign.setSeq1(comp_sequence);
  std::cout<<"start extensions "<<std::flush;
  for(MapPath::iterator iter_hash=map_path.begin();iter_hash!=map_path.end();
      iter_hash++)
    {
      SDGBioSeq sseq=subject_db[iter_hash->first.second-1];
      fastExtAlign.setSeq2(sseq);
      rfastExtAlign.setSeq2(sseq);
      for(std::list<RangePairSet>::iterator i=iter_hash->second.begin();
	  i!=iter_hash->second.end();i++)
		{
		  if(i->getRangeQ().getLength()>min_size
			 && i->getRangeS().getLength()>min_size)
			{
	 	      std::cout<<"."<<std::flush;
	 	     if(verbose>0)
	 	     {
	 	    	 std::cout<<"Query length:"<<sequence.length()<<std::endl;
	 	    	 std::cout<<"Subject length:"<<sseq.length()<<std::endl;
	 	    	 std::cout<<*i<<std::endl;
	 	     }
			  if(i->getRangeQ().isPlusStrand()
			 && i->getRangeS().isPlusStrand())
				{
				  fastExtAlign.setStart(i->getRangeQ().getEnd(),
							i->getRangeS().getEnd(),extend_len);
				  fastExtAlign.extend_dir(i->getScore());
				  if(verbose>0)
				  {
					  std::cout<<"extended end by "<<fastExtAlign.getEndSeq1()-i->getRangeQ().getEnd()<<std::endl;
				  }
				  i->getRangeQ().setEnd(fastExtAlign.getEndSeq1());
				  i->getRangeS().setEnd(fastExtAlign.getEndSeq2());

				  fastExtAlign.setStart(i->getRangeQ().getStart(),
							i->getRangeS().getStart(),extend_len);
				  fastExtAlign.extend_rev(i->getScore());
				  if(verbose>0)
				  {
					  std::cout<<"extended start by "<<i->getRangeQ().getStart()-fastExtAlign.getEndSeq1()<<std::endl;
				  }
				  i->getRangeQ().setStart(fastExtAlign.getEndSeq1());
				  i->getRangeS().setStart(fastExtAlign.getEndSeq2());
				}
			  if( !i->getRangeQ().isPlusStrand()
			  && i->getRangeS().isPlusStrand())
				{
				  unsigned qs=seqlen-i->getRangeQ().getStart()+1;
				  unsigned qe=seqlen-i->getRangeQ().getEnd()+1;
				  rfastExtAlign.setStart(qe,i->getRangeS().getEnd(),
							extend_len);
				  rfastExtAlign.extend_dir(i->getScore());
				  if(verbose>0)
				  {
					  std::cout<<"extended end by "<<rfastExtAlign.getEndSeq1()-i->getRangeQ().getEnd()<<std::endl;
				  }
				  i->getRangeQ().setEnd(seqlen
							-rfastExtAlign.getEndSeq1()+1);
				  i->getRangeS().setEnd(rfastExtAlign.getEndSeq2());

				  rfastExtAlign.setStart(qs,i->getRangeS().getStart(),
							extend_len);
				  rfastExtAlign.extend_rev(i->getScore());
				  if(verbose>0)
				  {
					  std::cout<<"extended start by "<<i->getRangeQ().getStart()-rfastExtAlign.getEndSeq1()<<std::endl;
				  }
				  i->getRangeQ().setStart(seqlen
							  -rfastExtAlign.getEndSeq1()+1);
				  i->getRangeS().setStart(rfastExtAlign.getEndSeq2());

				}
			}
		}
    }
  std::cout<<"ok"<<std::endl;
}





