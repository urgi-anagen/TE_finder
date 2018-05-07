#include "Hasher.h"

//-------------------------------------------------------------------------
void Hasher::load(const SDGString& filenameS)
{
  
  subject_db.load(filenameS);
  nbseqS=0;
  nb_word=0;
  if(!read_idx(filenameS))
    {
      std::cout<<"Counting words ... "<<std::flush;
      SDGFastaIstream inS(filenameS);
      for(unsigned k=0;k<max_key;k++)
    	  word_count[k]=0;
      
      while(inS)
		{
		  SDGBioSeq sS;
		  if(inS)
			inS>>sS;
		  //std::cout<<sS.getDE()<<std::endl;
		  hashSeqCount(sS);
		}
      inS.close();
      std::cout<<nb_word<<std::endl;

      std::cout<<"Stat:"<<std::endl;
      std::vector<unsigned> distr=word_count;
      sort(distr.begin(),distr.end());
      std::vector<unsigned> clean_distr;
      copy(lower_bound(distr.begin(),distr.end(),(unsigned)1),
	   distr.end(),back_inserter(clean_distr));
      std::cout<<" min="<<clean_distr.front()<<std::endl;
      std::cout<<" 1st quartile="
	       << clean_distr[(unsigned)floor(0.25*clean_distr.size())]<<std::endl;
      std::cout<<" median="
	       << clean_distr[(unsigned)floor(0.5*clean_distr.size())]<<std::endl;
      std::cout<<" 3st quartile="
	       << clean_distr[(unsigned)floor(0.75*clean_distr.size())]<<std::endl;
      std::cout<<" max="<< clean_distr.back()<<std::endl;
      unsigned max_count=
	clean_distr[clean_distr.size()-(unsigned)(cutoff*clean_distr.size())-1];
      std::cout<<"=>cut-off="<<cutoff<<" words occuring more than "<<max_count
	       <<" are removed!"<<std::endl;
      for(unsigned k=0;k<max_key;k++)
    	  if(word_count[k]>max_count)
    		  word_count[k]=0;
      word_count[0]=0; //remove words AAAAAA... NNNNNN.... XXXXX....
      std::cout<<"Load word positions"<<std::endl;
      word_pos.resize(nb_word);
      
      std::cout<<"Prepare hash table pointers"<<std::endl;
      unsigned k=0;
      std::vector<WordSpos>::iterator last_it=word_pos.begin();
      for(std::vector<unsigned>::iterator i=word_count.begin(); i!=word_count.end();i++)
		{
		  hash2wpos[k]=last_it;
		  last_it=last_it+(*i);
		  k++;
		}
      hash2wpos[k]=last_it;

      hash_ptr=hash2wpos;
      SDGFastaIstream inS2(filenameS);
      while(inS2)
		{
		  SDGBioSeq sS;
		  if(inS2)
			inS2>>sS;
		  hashSeqPos(sS);
		}
      inS2.close();
      hash_ptr.clear();
      save_idx(filenameS);
    }
  std::cout<<"end preparation of subjects"<<std::endl;
}
//-------------------------------------------------------------------------
// Read a hash index
bool Hasher::read_idx(const SDGString& filename)
{
  std::ifstream in(filename+".idx");
  if(!in) return false;
  std::cout<<"Read idx"<<std::endl;
  word_pos.clear();
  unsigned val;
  in>>val;
  if(val!=word_size) return false;
  in>>val;
  for(unsigned i=0;i<val;i++)
    {
      unsigned p, n;
      in>>p>>n;
      WordSpos wpos(p,n);
      word_pos.push_back(wpos);
    }
  word_count.clear();
  in>>val;
  for(unsigned i=0;i<val;i++)
    {
      unsigned n;
      in>>n;
      word_count.push_back(n);
    }
  in.close();

  std::cout<<"Prepare hash table pointers"<<std::endl;
  unsigned k=0;
  std::vector<WordSpos>::iterator last_it=word_pos.begin();
  for(std::vector<unsigned>::iterator i=word_count.begin(); i!=word_count.end();i++)
    {
      hash2wpos[k]=last_it;
      last_it=last_it+(*i);
      k++;
    }
  hash2wpos[k]=last_it;

  return true;
}
//-------------------------------------------------------------------------
// Save a hash index
void Hasher::save_idx(const SDGString& filename)
{
  std::ofstream fout(filename+".idx");
  std::cout<<"Save idx"<<std::endl;
  fout<<word_size<<std::endl;
  fout<<word_pos.size()<<std::endl;
  for(std::vector<WordSpos>::iterator i=word_pos.begin();i!=word_pos.end();i++)
    {
      fout<<i->pos<<" "<<i->numSeq<<"\t";
    }
  fout<<word_count.size()<<std::endl;
  for(std::vector<unsigned>::iterator i=word_count.begin();
      i!=word_count.end();i++)
    {
      fout<<*i<<"\t";
    }
}
//-------------------------------------------------------------------------
// Count words
void Hasher::hashSeqCount(const SDGBioSeq& seq)
{
  unsigned len=seq.length();
  std::string::const_iterator s=seq.begin();
  unsigned last_pos=len-word_size;
  for(unsigned i=0;i<=last_pos;i+=step_s)
    {
      word_count[hseq(s)]++;
      nb_word++;
      s+=step_s;
    }
}
//-------------------------------------------------------------------------
void Hasher::hashSeqPos(const SDGBioSeq& seq)
{
  nbseqS++;
  unsigned len=seq.length();
  std::string::const_iterator s=seq.begin();
  unsigned last_pos=len-word_size;
  unsigned key;
  for(unsigned i=0;i<=last_pos;i+=step_s)
    {
      key=hseq(s);

      if(word_count[key]!=0)
		{
		  *(hash_ptr[key])=WordSpos(i,nbseqS);
		  hash_ptr[key]++;
		}
      s+=step_s;
    }
}
//-------------------------------------------------------------------------
void Hasher::search(const SDGBioSeq& seq)
{
  std::cout<<"start hashing..."<<std::flush;
  nbseqQ++;
  sequence=seq;
  comp_sequence=seq.complement();
  if(repeat_mode) matchWords_repeat();
  else matchWords();
  std::cout<<"ok"<<std::endl;
  std::cout<<"start word connexion..."<<std::flush;
  diagSearch();
  std::cout<<"ok"<<std::endl;

}
//-------------------------------------------------------------------------
// Search for alignments with word matches
void Hasher::matchWords(void)
{
  diag_map.clear();
  diag_map_comp.clear();
  unsigned last_pos=sequence.length()-step_q;

  std::string::const_iterator seq=sequence.begin();
  std::string::const_iterator seq_comp=comp_sequence.begin();

  unsigned key_d=0,key_c=0;
  for(unsigned i=0;i<=last_pos;i+=step_q)
    {
      key_d=hseq(seq);
      seq+=step_q;

      std::vector<WordSpos>::iterator begin_d=hash2wpos[key_d];
      std::vector<WordSpos>::iterator end_d=hash2wpos[key_d+1];
      for(std::vector<WordSpos>::iterator j=begin_d;j!=end_d;j++)
	  diag_map.push_back(Diag(i-j->pos,j->pos,j->numSeq));

      key_c=hseq(seq_comp);
      seq_comp+=step_q;

      std::vector<WordSpos>::iterator begin_c=hash2wpos[key_c];
      std::vector<WordSpos>::iterator end_c=hash2wpos[key_c+1];
      for(std::vector<WordSpos>::iterator j=begin_c;j!=end_c;j++)
	  diag_map_comp.push_back(Diag(i-j->pos,j->pos,j->numSeq));

    }
}
//-------------------------------------------------------------------------
// Search for repeats with word matches
void Hasher::matchWords_repeat(void)
{
  diag_map.clear();
  diag_map_comp.clear();
  unsigned last_pos=sequence.length()-step_q;

  std::string::const_iterator seq=sequence.begin();
  std::string::const_iterator seq_comp=comp_sequence.begin();

  unsigned key_d=0,key_c=0;
  for(unsigned i=0;i<=last_pos;i+=step_q)
    {
      key_d=hseq(seq);
      seq+=step_q;

      std::vector<WordSpos>::iterator begin_d=hash2wpos[key_d];
      std::vector<WordSpos>::iterator end_d=hash2wpos[key_d+1];
      for(std::vector<WordSpos>::iterator j=begin_d;j!=end_d;j++)
      {
    	  if(j->numSeq==nbseqQ && i==j->pos)
    		  continue;
    	  diag_map.push_back(Diag(i-j->pos,j->pos,j->numSeq));
      }

      key_c=hseq(seq_comp);
      seq_comp+=step_q;

      std::vector<WordSpos>::iterator begin_c=hash2wpos[key_c];
      std::vector<WordSpos>::iterator end_c=hash2wpos[key_c+1];
      for(std::vector<WordSpos>::iterator j=begin_c;j!=end_c;j++)
          for(std::vector<WordSpos>::iterator j=begin_d;j!=end_d;j++)
          {
//        	  if(j->numSeq==nbseqQ && i==j->pos)
//        		  continue;
        	  diag_map_comp.push_back(Diag(i-j->pos,j->pos,j->numSeq));
          }
    }
}
//-------------------------------------------------------------------------
// Search for diagonal of word matches
void Hasher::diagSearch(void)
{
  frag.clear();
  unsigned count=0;
  unsigned size=diag_map.size(); 
  unsigned len=sequence.length();
  unsigned wdlen=wdist*word_size;
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
	      RangeAlign r1(nbseqQ,diag+start+1,diag+end+word_size);
	      RangeAlign r2(numSeq,start+1,end+word_size);
	      
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
	  RangeAlign r1(nbseqQ,diag+start+1,diag+end+word_size);
	  RangeAlign r2(numSeq,start+1,end+word_size);
	      
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
			    len-(diag+end)-word_size+1);
	      RangeAlign r2(numSeq,start+1,end+word_size);
	      
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
			    len-(diag+end)-word_size+1);
	  RangeAlign r2(numSeq,start+1,end+word_size);
	  
	  
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
void Hasher::fragAlign(double match,double mism, double gopen,
			   double gext, unsigned over, bool join)
{
  std::cout<<"start fragment connexion ..."<<std::flush;
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
void Hasher::print_frag(std::ostream& out)
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
void Hasher::print(unsigned min_size,std::ostream& out)
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
void Hasher::write(unsigned min_size,std::ostream& out)
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
void Hasher::write_align(unsigned min_size,std::ostream& out)
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
void Hasher::extend(unsigned min_size)
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
			  if(i->getRangeQ().isPlusStrand()
			 && i->getRangeS().isPlusStrand())
				{
				  fastExtAlign.setStart(i->getRangeQ().getEnd(),
							i->getRangeS().getEnd(),extend_len);
				  fastExtAlign.extend_dir(i->getScore());
				  //std::cout<<"extended end by "<<fastExtAlign.getEndSeq1()-i->getRangeQ().getEnd()<<std::endl;
				  i->getRangeQ().setEnd(fastExtAlign.getEndSeq1());
				  i->getRangeS().setEnd(fastExtAlign.getEndSeq2());

				  fastExtAlign.setStart(i->getRangeQ().getStart(),
							i->getRangeS().getStart(),extend_len);
				  fastExtAlign.extend_rev(i->getScore());
				  //std::cout<<"extended start by "<<i->getRangeQ().getStart()-fastExtAlign.getEndSeq1()<<std::endl;
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
				  i->getRangeQ().setEnd(seqlen
							-rfastExtAlign.getEndSeq1()+1);
				  i->getRangeS().setEnd(rfastExtAlign.getEndSeq2());

				  rfastExtAlign.setStart(qs,i->getRangeS().getStart(),
							extend_len);
				  rfastExtAlign.extend_rev(i->getScore());
				  i->getRangeQ().setStart(seqlen
							  -rfastExtAlign.getEndSeq1()+1);
				  i->getRangeS().setStart(rfastExtAlign.getEndSeq2());

				}
			}
		}
    }
  std::cout<<"ok"<<std::endl;
}







