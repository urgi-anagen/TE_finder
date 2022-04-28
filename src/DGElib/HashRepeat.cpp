#include "HashRepeat.h"
//-------------------------------------------------------------------------
//void HashRepeat::search_reputer(SDGBioSeq& seq, unsigned errors)
//{
//  nbseqQ++;
//  queryName.push_back(seq.getDE());
//  sequence=seq;
//  reputer(errors);
//  fragAlign(1,0.8,1.6,0.4,0);
//  extend();
//}
////-------------------------------------------------------------------------
//void HashRepeat::reputer(unsigned errors)
//{
//  SDGString tmp_filename=SDGString(getpid())+".tmp";
//  SDGFastaOstream out(tmp_filename);
//  out<<sequence;
//  out.close();
//
//  std::ostringstream cmd;
//  if(errors==0)
//    cmd<<"repfind -f -p -l "<<word_size
//       <<" "<<tmp_filename<<" > "<<tmp_filename<<".res";
//  else
//    cmd<<"repfind -f -p -l "<<word_size
//       <<" -e "<<errors
//       <<" "<<tmp_filename<<" > "<<tmp_filename<<".res";
//
//  system(cmd.str().c_str());
//
//  std::ostringstream filename_res;
//  filename_res<<tmp_filename<<".res";
//
//  frag.clear();
//  std::ifstream in(filename_res.str().c_str());
//  unsigned count=0;
//  while(in)
//    {
//      RangePair rp;
//      rp.readReputer(in);
//      if(rp.getLength())
//	{
//	  rp.setId(++count);
//	  rp_list.push_back(rp);
//	}
//    }
//  in.close();
//
//  std::ostringstream rm_cmd;
//  rm_cmd<<"rm "<<tmp_filename<<"*";
//  system(rm_cmd.str().c_str());
//}
//-------------------------------------------------------------------------
void HashRepeat::search(SDGBioSeq& seq)
{
  nbseqQ++;
  queryName.push_back(seq.getDE());
  sequence=seq;
  sequence_comp=sequence.complement();
  prepare();
  matchWords();
  diagSearch();
  fragAlign(1,0.8,1.6,0.4,0);
  extend();
}
//-------------------------------------------------------------------------
void HashRepeat::prepare(double cutoff)
{
  nb_word=0;
  std::string::const_iterator s=sequence.begin();
  unsigned last_pos=sequence.length()-word_size;
  fill(word_count.begin(),word_count.end(),0);
  for(unsigned i=0;i<=last_pos;i+=word_size)
    {
      word_count[hseq(s)]++;
      nb_word++;
      s+=word_size;
    }
  std::cout<<"Stat:";
  std::vector<unsigned> distr=word_count;
  sort(distr.begin(),distr.end());
  std::vector<unsigned> clean_distr;
  copy(lower_bound(distr.begin(),distr.end(),(unsigned)1),
       distr.end(),back_inserter(clean_distr));
  std::cout<<" min="<<clean_distr.front();
  std::cout<<" 1st quartile="
      << clean_distr[(unsigned)floor(0.25*clean_distr.size())];
  std::cout<<" median="
      << clean_distr[(unsigned)floor(0.5*clean_distr.size())];
  std::cout<<" 3st quartile="
      << clean_distr[(unsigned)floor(0.75*clean_distr.size())];
  std::cout<<" max="<< clean_distr.back()<<std::endl;
  unsigned max_count=
    clean_distr[clean_distr.size()-(unsigned)(cutoff*clean_distr.size())-1];
  std::cout<<"# words="<<clean_distr.size()<<"=>cut-off="<<cutoff<<" words occuring more than "<<max_count
      <<" are removed!"<<std::endl;
  for(unsigned k=0;k<max_key;k++)
    if(word_count[k]>max_count)
      word_count[k]=0;
  word_count[0]=0; //remove_self_hits words AAAAAA... NNNNNN.... XXXXX....
  word_pos.resize(nb_word);
  unsigned k=0;
  std::vector<unsigned>::iterator last_it=word_pos.begin();
  for(std::vector<unsigned>::iterator i=word_count.begin();
      i!=word_count.end();i++)
    {
      hash2wpos[k]=last_it;
      last_it=hash2wpos[k]+(*i);
      k++;
    }
  hash2wpos[k]=word_pos.end();
  hash_ptr=hash2wpos;
}
//-------------------------------------------------------------------------
void HashRepeat::matchWords(void)
{
  diag_map.clear();
  diag_map_comp.clear();
  unsigned len=sequence.length();
  unsigned last_pos=len-word_size;

  std::string::const_iterator seq=sequence.begin();

  unsigned key_d=0;
  for(unsigned i=0;i<=last_pos;i++)
    {
      key_d=hseq(seq);
      seq++;

      std::vector<unsigned>::iterator begin_d=hash2wpos[key_d];
      std::vector<unsigned>::iterator end_d=hash_ptr[key_d];
      for(std::vector<unsigned>::iterator j=begin_d;j!=end_d;j++)
	diag_map.push_back(Diag(i-(*j),(*j)));

      if(i%word_size==0 && word_count[key_d]!=0)
	{
	  *(hash_ptr[key_d])=i;
	  hash_ptr[key_d]++;
	}      
    }

  std::string::const_iterator seq_comp=sequence_comp.begin();
  unsigned key_c=0;
  for(unsigned i=0;i<=last_pos;i++)
    {
      key_c=hseq(seq_comp);
      seq_comp++;

      std::vector<unsigned>::iterator begin_c=hash2wpos[key_c];
      std::vector<unsigned>::iterator end_c=hash_ptr[key_c];
      unsigned rev_i=len-i;
      for(std::vector<unsigned>::iterator j=begin_c;j!=end_c && *j<rev_i;j++)
	diag_map_comp.push_back(Diag(i-(*j),(*j)));

    }
}
//-------------------------------------------------------------------------
void HashRepeat::diagSearch(void)
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
      int diag=0;
      bool in=false; 
      bool found=false;

      Diag& prev_d=diag_map[0];
      for( unsigned i=1; i<size; ++i)
	{
	  Diag& curr_d=diag_map[i];
	  if(prev_d.diag==curr_d.diag)
	    {
	      if(prev_d.wpos+word_size>=curr_d.wpos)
		{
		  if(in)
		    {
		      end=curr_d.wpos;
		    }
		  else
		    {
		      diag=prev_d.diag;
		      start=prev_d.wpos;
		      end=curr_d.wpos;
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
	      RangeAlign r2(nbseqQ,start+1,end+word_size);
	      
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
	  RangeAlign r2(nbseqQ,start+1,end+word_size);
	      
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
      int diag=0;
      bool in=false; 
      bool found=false;
      Diag& prev_d=diag_map_comp[0];
      for( unsigned i=1; i<size; ++i)
	{
	  Diag& curr_d=diag_map_comp[i];
	  if(prev_d.diag==curr_d.diag) 
	    {
	      if(prev_d.wpos+word_size>=curr_d.wpos)
		{
		  if(in)
		    {
		      end=curr_d.wpos;
		    }
		  else
		    {
		      diag=prev_d.diag;
		      start=prev_d.wpos;
		      end=curr_d.wpos;
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
	      RangeAlign r2(nbseqQ,start+1,end+word_size);
	      
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
	  RangeAlign r2(nbseqQ,start+1,end+word_size);
	  
	  
	  RangePair rp(r1,r2);
	  rp.setId(++count);
	  rp.setScore(r1.getLength());
	  rp.setIdentity(1.00);      
	  rp.setLength(r1.getLength());
	  frag.push_back(std::pair<unsigned,RangePair>(rp.getScore(),rp));	  
	}
    }

  rp_list.clear();

  for(std::vector< std::pair<unsigned,RangePair> >::iterator i=frag.begin();
      i!=frag.end();i++)
    {
      if(i->second.getRangeQ().getLength()>=2*word_size)
	rp_list.push_back(i->second);
    }
}
//-------------------------------------------------------------------------
void HashRepeat::fragAlign(double match,double mism, double gopen,
			   double gext, unsigned over)
{
  FragAlign fragAligner(mism,gopen,gext,over);

  path=fragAligner.join(rp_list);
}
//-------------------------------------------------------------------------
void HashRepeat::extend(void)
{
  FastExtAlign fastExtAlign;
  FastExtAlign rfastExtAlign;
  unsigned seqlen=sequence.length();

  fastExtAlign.setSeq1(sequence);
  rfastExtAlign.setSeq1(sequence_comp);
  fastExtAlign.setSeq2(sequence);
  rfastExtAlign.setSeq2(sequence);
  for(std::list<RangePairSet>::iterator i=path.begin();i!=path.end();
      i++)
    {
      if(i->getRangeQ().isPlusStrand() 
	 && i->getRangeS().isPlusStrand())
	{
	  fastExtAlign.setStart(i->getRangeQ().getEnd(),
				i->getRangeS().getEnd(),extend_len);
	  int score1=fastExtAlign.extend_dir(i->getScore());	  
	  i->getRangeQ().setEnd(fastExtAlign.getEndSeq1());
	  i->getRangeS().setEnd(fastExtAlign.getEndSeq2());
	      
	  fastExtAlign.setStart(i->getRangeQ().getStart(),
				i->getRangeS().getStart(),extend_len);
	  int score2=fastExtAlign.extend_rev(i->getScore());

	  i->getRangeQ().setStart(fastExtAlign.getEndSeq1());
	  i->getRangeS().setStart(fastExtAlign.getEndSeq2());	      
	  i->setScore(score1+score2-(i->getScore()));
	}
      if( !i->getRangeQ().isPlusStrand() 
	  && i->getRangeS().isPlusStrand())
	{
	  unsigned qs=seqlen-i->getRangeQ().getStart()+1;
	  unsigned qe=seqlen-i->getRangeQ().getEnd()+1;
	  rfastExtAlign.setStart(qe,i->getRangeS().getEnd(),
				 extend_len);
	  int score1=rfastExtAlign.extend_dir(i->getScore());	  
	  i->getRangeQ().setEnd(seqlen
				-rfastExtAlign.getEndSeq1()+1);
	  i->getRangeS().setEnd(rfastExtAlign.getEndSeq2());
	  
	  rfastExtAlign.setStart(qs,i->getRangeS().getStart(),
				 extend_len);
	  int score2=rfastExtAlign.extend_rev(i->getScore());	  
	  i->getRangeQ().setStart(seqlen
				  -rfastExtAlign.getEndSeq1()+1);
	  i->getRangeS().setStart(rfastExtAlign.getEndSeq2());
	  
	  i->setScore(score1+score2-(i->getScore()));
	}
    }
}
//-------------------------------------------------------------------------
void HashRepeat::print(std::ostream& out)
{
  out<<"\nFragments found:"<<std::endl;
  for(std::list<RangePair>::iterator i=rp_list.begin();
      i!=rp_list.end();i++)
    {
      out<<queryName[i->getRangeQ().getNumChr()-1]<<"\t"
	 <<i->getRangeQ().getStart()<<".."<<i->getRangeQ().getEnd()
	 <<"\t"
	 <<queryName[i->getRangeS().getNumChr()-1]<<"\t"
	 <<i->getRangeS().getStart()<<".."<<i->getRangeS().getEnd()
	 <<"\t"
	 <<i->getScore()<<std::endl;
    } 
}
//-------------------------------------------------------------------------
void HashRepeat::write(std::ostream& out)
{
  for(std::list<RangePairSet>::iterator i=path.begin();
      i!=path.end();i++)
    {
      out<<queryName[i->getRangeQ().getNumChr()-1]<<"\t"
	 <<i->getRangeQ().getStart()<<"\t"<<i->getRangeQ().getEnd()
	 <<"\t"
	 <<queryName[i->getRangeS().getNumChr()-1]<<"\t"
	 <<i->getRangeS().getStart()<<"\t"<<i->getRangeS().getEnd()
	 <<"\t"
	 <<i->getScore()<<std::endl;
    } 
}
//-------------------------------------------------------------------------
void HashRepeat::writeMap(std::ostream& out)
{
  unsigned count=0;
  for(std::list<RangePairSet>::iterator i=path.begin();
      i!=path.end();i++)
    {
      SDGString name,info;
      if(i->getRangeS().isPlusStrand())
	{
	  name="Rep"+SDGString(count++)+" Dir";
	}
      else
	{
	  name="Rep"+SDGString(count++)+" Rev";
	}	
      
      info=" | size="+SDGString(i->getRangeQ().getLength());
      info+=" | score="+SDGString(i->getScore());
      info+=" | lenTE="+SDGString(i->getRangeS().getMax()
				  -i->getRangeQ().getMin()+1);
      out<<name+info<<"\t"<<sequence.getDE()<<"\t"
	 <<i->getRangeQ().getStart()<<"\t"<<i->getRangeQ().getEnd()<<std::endl;
      
      info=" | size="+SDGString(i->getRangeS().getLength());
      info+=" | score="+SDGString(i->getScore());
      info+=" | lenTE="+SDGString(i->getRangeS().getMax()
				  -i->getRangeQ().getMin()+1);
      out<<name+info<<"\t"<<sequence.getDE()<<"\t"
	 <<i->getRangeS().getStart()<<"\t"<<i->getRangeS().getEnd()<<std::endl;
    } 
}
//-------------------------------------------------------------------------
unsigned HashRepeat::writeSet(std::ostream& out,unsigned count)
{
  for(std::list<RangePairSet>::iterator i=path.begin();
      i!=path.end();i++)
    {
      SDGString name,info;
      unsigned id=++count;
      if(i->getRangeS().isPlusStrand())
	name="RepDir";
      else
	name="RepRev";
      
      
      info="|size="+SDGString(i->getRangeQ().getLength());
      info+="|score="+SDGString(i->getScore());
      info+="|lenTE="+SDGString(i->getRangeS().getMax()
				-i->getRangeQ().getMin()+1);
      
      if(i->getRangeS().isPlusStrand())
	out<<id<<"\t"<<name+info<<"\t"<<sequence.getDE()<<"\t"
	   <<i->getRangeQ().getMin()<<"\t"
	   <<i->getRangeQ().getMax()<<std::endl;
      else
	out<<id<<"\t"<<name+info<<"\t"<<sequence.getDE()<<"\t"
	   <<i->getRangeQ().getMax()<<"\t"
	   <<i->getRangeQ().getMin()<<std::endl;
      
      info="|size="+SDGString(i->getRangeS().getLength());
      info+="|score="+SDGString(i->getScore());
      info+="|lenTE="+SDGString(i->getRangeS().getMax()
				-i->getRangeQ().getMin()+1);
      if(i->getRangeS().isPlusStrand())
	out<<id<<"\t"<<name+info<<"\t"<<sequence.getDE()<<"\t"
	   <<i->getRangeS().getMin()<<"\t"
	   <<i->getRangeS().getMax()<<std::endl;
      else
	out<<id<<"\t"<<name+info<<"\t"<<sequence.getDE()<<"\t"
	   <<i->getRangeS().getMax()<<"\t"
	   <<i->getRangeS().getMin()<<std::endl;
    }
  return count;
}
// //-------------------------------------------------------------------------
// void HashRepeat::plot(const SDGString& title)
// {
//   if(!h) return;
//   std::ostringstream cmd;
//   cmd<<"set title \" "<<title<<" \""<<std::endl
//      <<"set xrange [1:"<<sequence.length()<<"]"<<std::endl
//      <<"set yrange [1:"<<sequence.length()<<"]"<<std::endl;

//   bool frag_dir=false,frag_rev=false,join=false;
//   unsigned max_x=0,max_y=0;
//   unsigned min_x=UINT_MAX,min_y=UINT_MAX;

//   std::ofstream fout1(datafile1);
//   std::ofstream fout2(datafile2);
//   std::ofstream fout3(datafile3);

//   for(std::list<RangePairSet>::iterator k=path.begin();
//       k!=path.end();k++)
//     {
//       for(std::list<RangePair>::iterator j=k->begin();
// 	  j!=k->end();j++)
// 	{
// 	  if(max_x<
// 	     std::max(j->getRangeQ().getStart(),j->getRangeQ().getEnd()) )
// 	    max_x=std::max(j->getRangeQ().getStart(),
// 		      j->getRangeQ().getEnd());
// 	  if(max_y<
// 	     std::max(j->getRangeS().getStart(),j->getRangeS().getEnd()) )
// 	    max_y=std::max(j->getRangeS().getStart(),
// 		      j->getRangeS().getEnd());
// 	  if(min_x>
// 	     std::min(j->getRangeQ().getStart(),j->getRangeQ().getEnd()) )
// 	    min_x=std::min(j->getRangeQ().getStart(),
// 		      j->getRangeQ().getEnd());
// 	  if(min_y>
// 	     std::min(j->getRangeS().getStart(),j->getRangeS().getEnd()) )
// 	    min_y=std::min(j->getRangeS().getStart(),
// 		      j->getRangeS().getEnd());


// 	  if(j->getRangeQ().isPlusStrand() 
// 	     && j->getRangeS().isPlusStrand())
// 	    { 
// 	      fout1<<j->getRangeQ().getStart()
// 		   <<"\t"<<j->getRangeS().getStart()
// 		   <<std::endl;
// 	      fout1<<j->getRangeQ().getEnd()
// 		   <<"\t"<<j->getRangeS().getEnd()
// 		   <<std::endl;
// 	      fout1<<std::endl;
// 	      frag_dir=true;
// 	    }
// 	  else
// 	    {
// 	      fout2<<j->getRangeQ().getStart()
// 		   <<"\t"<<j->getRangeS().getStart()
// 		   <<std::endl;
// 	      fout2<<j->getRangeQ().getEnd()
// 		   <<"\t"<<j->getRangeS().getEnd()
// 		   <<std::endl;
// 	      fout2<<std::endl;
// 	      frag_rev=true;
// 	    }
// 	}
      
//       for(std::list<RangePair>::iterator j=k->begin();
// 	  j!=k->end();j++)
// 	{
// 	  std::list<RangePair>::iterator l=j;
// 	  if(++l==k->end()) break;
// 	  if(j->getRangeQ().isPlusStrand())
// 	    {
// 	      fout3<<std::max(j->getRangeQ().getEnd(),
// 			 j->getRangeQ().getStart());
// 	    }
// 	  else
// 	    {
// 	      fout3<<std::min(j->getRangeQ().getEnd(),
// 			 j->getRangeQ().getStart());
// 	    }
// 	  if(j->getRangeS().isPlusStrand())
// 	    {
// 	      fout3<<"\t"<<std::max(j->getRangeS().getEnd(),
// 			       j->getRangeS().getStart())<<std::endl;
// 	    }
// 	  else
// 	    {
// 	      fout3<<"\t"<<std::min(j->getRangeS().getEnd(),
// 			       j->getRangeS().getStart())<<std::endl;
// 	    }
// 	  if(j->getRangeQ().isPlusStrand())
// 	    {
// 	      fout3<<std::min(l->getRangeQ().getStart(),
// 			 l->getRangeQ().getEnd());
// 	    }
// 	  else
// 	    {
// 	      fout3<<std::max(l->getRangeQ().getStart(),
// 			 l->getRangeQ().getEnd());
// 	    }
// 	  if(j->getRangeS().isPlusStrand())
// 	    {
// 	      fout3<<"\t"<<std::min(l->getRangeS().getStart(),
// 			       l->getRangeS().getEnd())<<std::endl;
// 	    }
// 	  else
// 	    {
// 	      fout3<<"\t"<<std::max(l->getRangeS().getStart(),
// 			       l->getRangeS().getEnd())<<std::endl;
// 	    }
// 	  fout3<<std::endl;
// 	  join=true;
// 	}
//     }
//   fout1.close();
//   fout2.close();
//   fout3.close();
//   bool first=true;
//   if(path.size())
//     {
//       cmd<<" plot \"";
//       if(frag_dir)
// 	{
// 	  cmd<<datafile1<<"\" notitle with lines lt 1";
// 	  first=false;
// 	}
//       if(frag_rev)
// 	{
// 	  if(!first)
// 	    cmd<<", \"";
// 	  cmd<<datafile2<<"\" notitle with lines lt 1";
// 	  first=false;
// 	}
//       if(join)
// 	cmd<<", \""<<datafile3<<"\" notitle with lines lt 3";
//     }
//   gnuplot_cmd(h,(char*)cmd.str().c_str()) ;
// }
//-------------------------------------------------------------------------
//double HashRepeat::compress(unsigned min_size, double cover)
//{
//
//  RangeGroups grouper(std::numeric_limits<unsigned int>::max(),cover);
//
//  std::list<RangePair> signif_rp;
//  for(std::list<RangePairSet>::iterator i=path.begin();
//      i!=path.end();i++)
//    if(i->getRangeQ().getLength()>= min_size
//       && i->getRangeS().getLength()>= min_size)
//      signif_rp.push_back(RangePair(*i));
//
//  grouper.grouping(signif_rp);
//
//  GROUPLIST group_list=grouper.getGroupList();
//
//  unsigned sum=0,gr=0;
//  for(GROUPLIST::iterator i=group_list.begin();
//      i!=group_list.end();i++)
//    {
//      gr++;
//      for(LISTMEMBER::iterator j=i->begin();
//	  j!=i->end();j++)
//	{
//	  sum+=j->getLength();
//	}
//    }
//  return (double)sum/sequence.length();
//}






