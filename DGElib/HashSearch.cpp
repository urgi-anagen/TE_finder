#include "HashSearch.h"

//-------------------------------------------------------------------------
void HashSearch::load(SDGString& filenameS)
{
  
  subjectName.clear();
  nbseqS=0;
  nb_word=0;
  std::cout<<"Counting words ... "<<std::flush;
  SDGFastaIstream inS(filenameS);
  while(inS)
    {
      SDGBioSeq sS;
      if(inS)
	inS>>sS;
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
  std::cout<<"Prepare hash table pointers"<<std::endl;
  word_pos.resize(nb_word);
  unsigned k=0;
  std::vector<WordSpos>::iterator last_it=word_pos.begin();
  for(std::vector<unsigned>::iterator i=word_count.begin(); i!=word_count.end();i++)
    {
      hash2wpos[k]=last_it;
      last_it=hash2wpos[k]+(*i);
      k++;
    }
  hash2wpos[k]=word_pos.end();

  std::cout<<"Load word positions"<<std::endl;
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
  std::cout<<"end preparation of subjects"<<std::endl;

}
//-------------------------------------------------------------------------
void HashSearch::hashSeqCount(SDGBioSeq& seq)
{
  subjectName.push_back(seq.getDE());

  unsigned len=seq.length();
  std::string::const_iterator s=seq.begin();
  unsigned last_pos=len-word_size;
  for(unsigned i=0;i<=last_pos;i+=word_size)
    {
      word_count[hseq(s,word_size)]++;
      nb_word++;
      s+=word_size;
    }
}
//-------------------------------------------------------------------------
void HashSearch::hashSeqPos(SDGBioSeq& seq)
{
  nbseqS++;
  unsigned len=seq.length();
  std::string::const_iterator s=seq.begin();
  unsigned last_pos=len-word_size;
  unsigned key;
  for(unsigned i=0;i<=last_pos;i+=word_size)
    {
      key=hseq(s,word_size);;
      if(word_count[key]!=0)
	{
	  *(hash_ptr[key])=WordSpos(i+1,nbseqS);
	  hash_ptr[key]++;
	}

      s+=word_size;
    }
}
//-------------------------------------------------------------------------
void HashSearch::search(SDGBioSeq& seq)
{
  nbseqQ++;
  queryName.push_back(seq.getDE());
  sequence=seq;
  if(repeat_mode) matchWords_repeat();
  else matchWords();
  diagSearch();
}
//-------------------------------------------------------------------------
void HashSearch::matchWords(void)
{
  diag_map.clear();
  diag_map_comp.clear();
  unsigned last_pos=sequence.length()-word_size;

  std::string::const_iterator seq=sequence.begin();
  SDGBioSeq sequence_comp=sequence.complement();
  std::string::const_iterator seq_comp=sequence_comp.begin();

  unsigned key_d=0,key_c=0;
  for(unsigned i=1;i<=last_pos;i++)
    {
      if(i==1)
	{
	  key_d=hseq(seq,word_size);
	  seq+=word_size;
	}
      else
	key_d=hseq(seq++,1,key_d);

      std::vector<WordSpos>::iterator begin_d=hash2wpos[key_d];
      std::vector<WordSpos>::iterator end_d=hash2wpos[key_d+1];
      for(std::vector<WordSpos>::iterator j=begin_d;j!=end_d;j++)
	  diag_map.push_back(Diag(i-j->pos,j->pos,j->numSeq));

      if(i==1)
	{
	  key_c=hseq(seq_comp,word_size);
	  seq_comp+=word_size;
	}
      else
	key_c=hseq(seq_comp++,1,key_c);

      std::vector<WordSpos>::iterator begin_c=hash2wpos[key_c];
      std::vector<WordSpos>::iterator end_c=hash2wpos[key_c+1];
      for(std::vector<WordSpos>::iterator j=begin_c;j!=end_c;j++)
	  diag_map_comp.push_back(Diag(i-j->pos,j->pos,j->numSeq));

    }
}
//-------------------------------------------------------------------------
void HashSearch::matchWords_repeat(void)
{
  diag_map.clear();
  diag_map_comp.clear();
  unsigned last_pos=sequence.length()-word_size;

  std::string::const_iterator seq=sequence.begin();
  SDGBioSeq sequence_comp=sequence.complement();
  std::string::const_iterator seq_comp=sequence_comp.begin();

  unsigned key_d=0,key_c=0;
  for(unsigned i=1;i<=last_pos;i++)
    {
      if(i==1)
	{
	  key_d=hseq(seq,word_size);
	  seq+=word_size;
	}
      else
	key_d=hseq(seq++,1,key_d);

      std::vector<WordSpos>::iterator begin_d=hash2wpos[key_d];
      std::vector<WordSpos>::iterator end_d=hash2wpos[key_d+1];
      begin_d=lower_bound(begin_d,end_d,nbseqQ+1);
      for(std::vector<WordSpos>::iterator j=begin_d;j!=end_d;j++)
	    diag_map.push_back(Diag(i-j->pos,j->pos,j->numSeq));



      if(i==1)
	{
	  key_c=hseq(seq_comp,word_size);
	  seq_comp+=word_size;
	}
      else
	key_c=hseq(seq_comp++,1,key_c);
      std::vector<WordSpos>::iterator begin_c=hash2wpos[key_c];
      std::vector<WordSpos>::iterator end_c=hash2wpos[key_c+1];
      begin_c=lower_bound(begin_c,end_c,nbseqQ+1);
      for(std::vector<WordSpos>::iterator j=begin_c;j!=end_c;j++)
	diag_map_comp.push_back(Diag(i-j->pos,j->pos,j->numSeq));
    }
}
//-------------------------------------------------------------------------
void HashSearch::diagSearch(void)
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
	      if(prev_d.wpos.pos+word_size==curr_d.wpos.pos)
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
	      RangeAlign r1(numSeq,start,end+word_size);
	      RangeAlign r2(nbseqQ,diag+start,diag+end+word_size);
	      
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
	  RangeAlign r1(nbseqQ,start,end+word_size);
	  RangeAlign r2(numSeq,diag+start,diag+end+word_size);
	      
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
	      if(prev_d.wpos.pos+word_size==curr_d.wpos.pos)
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
	      RangeAlign r1(nbseqQ ,start,end+word_size-1);
	      
	      unsigned rstart,rend;
	      rstart=abs(diag-len+start);
	      rend=abs(diag-len+end+word_size-1);
	      
	      
	      RangeAlign r2(numSeq,rstart+1,rend+1);
	      
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
	  RangeAlign r1(nbseqQ ,start,end+word_size-1);
	  
	  unsigned rstart,rend;
	  rstart=abs(diag-len+start);
	  rend=abs(diag-len+end+word_size-1);
	  
	  
	  RangeAlign r2(numSeq,rstart+1,rend+1);
	  
	  RangePair rp(r1,r2);
	  rp.setId(++count);
	  rp.setScore(r1.getLength());
	  rp.setIdentity(1.00);      
	  rp.setLength(r1.getLength());

	  frag.push_back(std::pair<unsigned,RangePair>(rp.getScore(),rp));	  
	}
    }

  unsigned count2=0;
  rp_list.clear();

  sort(frag.begin(),frag.end(),comp());

  for(std::vector< std::pair<unsigned,RangePair> >::reverse_iterator i=frag.rbegin();
      i!=frag.rend() && count2<nfrag;i++)
    {
      count2++;
      rp_list.push_back(i->second);
    }
}
//-------------------------------------------------------------------------
void HashSearch::fragAlign(double match,double mism, double gopen,
			   double gext, unsigned over)
{
  FragAlign fragAligner(mism,gopen,gext,over);

  path=fragAligner.join(rp_list);
}
//-------------------------------------------------------------------------
void HashSearch::print(unsigned min_size,std::ostream& out, bool name)
{
  out<<"\nAlignments found:"<<std::endl;
  for(std::list<RangePair>::iterator i=rp_list.begin();
      i!=rp_list.end();i++)
    {
      if(i->getRangeQ().getLength()>min_size 
	 && i->getRangeS().getLength()>min_size)
	{
	  if(name)
	    {
	      out<<queryName[i->getRangeQ().getNumChr()-1]<<"\t"
		 <<i->getRangeQ().getStart()<<".."<<i->getRangeQ().getEnd()
		 <<"\t"
		 <<subjectName[i->getRangeS().getNumChr()-1]<<"\t"
		 <<i->getRangeS().getStart()<<".."<<i->getRangeS().getEnd()
		 <<"\t"
		 <<i->getScore()<<std::endl;
	    }
	  else
	    {
	      out<<i->getRangeQ().getNumChr()<<"\t"
		 <<i->getRangeQ().getStart()<<".."<<i->getRangeQ().getEnd()
		 <<"\t"
		 <<i->getRangeS().getNumChr()<<"\t"
		 <<i->getRangeS().getStart()<<".."<<i->getRangeS().getEnd()
		 <<"\t"
		 <<i->getScore()<<std::endl;
	    }
	}
    } 
}
//-------------------------------------------------------------------------
void HashSearch::write(unsigned min_size,std::ostream& out, bool name)
{
  for(std::list<RangePair>::iterator i=rp_list.begin();
      i!=rp_list.end();i++)
    {
      if(i->getRangeQ().getLength()>min_size 
	 && i->getRangeS().getLength()>min_size)
	{
	  if(name)
	    {
	      out<<queryName[i->getRangeQ().getNumChr()-1]<<"\t"
		 <<i->getRangeQ().getStart()<<"\t"<<i->getRangeQ().getEnd()
		 <<"\t"
		 <<subjectName[i->getRangeS().getNumChr()-1]<<"\t"
		 <<i->getRangeS().getStart()<<"\t"<<i->getRangeS().getEnd()
		 <<"\t"
		 <<i->getScore()<<std::endl;
	    }
	  else
	    {
	      out<<i->getRangeQ().getNumChr()<<"\t"
		 <<i->getRangeQ().getStart()<<"\t"<<i->getRangeQ().getEnd()
		 <<"\t"
		 <<i->getRangeS().getNumChr()<<"\t"
		 <<i->getRangeS().getStart()<<"\t"<<i->getRangeS().getEnd()
		 <<"\t"
		 <<i->getScore()<<std::endl;
	    }
	}
    } 
}







