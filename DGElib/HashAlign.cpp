#include "HashAlign.h"
//-------------------------------------------------------------------------
void HashAlign::load_matches(SDGString filename)
{
  std::ifstream in(filename);
  while(1)
    {
      RangePair rp;
      in>>rp;
      in.peek();
      if(!in)
	break;
      std::cout<<rp<<std::endl;
      sort_frag.insert(std::pair<unsigned long,RangePair>(rp.getScore(),rp));
    }
  in.close();
}
//-------------------------------------------------------------------------
void HashAlign::search(unsigned word_size)
{
  SDGString seq1=sequence1.toString();
  SDGString seq2=sequence2.toString();
  unsigned len1=seq1.length();
  unsigned len2=seq2.length();
  std::map<SDGString,std::list<unsigned> > map_word;
  std::map<int,Diag> diag_map;
  
  for(unsigned i=0;i<len1-word_size+1;i++)
    {
      SDGString word=seq1.substr(i,word_size);
      unsigned pos=word.rfind('N');
      if(pos>=word_size)
	  map_word[word].push_back(i);
    }
	
  for(unsigned i=0;i<len2-word_size+1;i++)
    {
      SDGString word=seq2.substr(i,word_size);
      unsigned pos=word.rfind('N');
      if(pos>=word_size)
	{
	  std::list<unsigned>& l=map_word[word];
	  for(std::list<unsigned>::iterator it=l.begin();
	      it!=l.end();it++)
	    {
	      diag_map[i-(*it)].add(Range((*it),(*it)+word_size-1));
	    }
	}
    }
  
  sort_frag.clear();
  
  for(std::map<int,Diag>::iterator i=diag_map.begin();
      i!=diag_map.end();i++)
    {
      for(Diag::iterator j=i->second.begin();
	  j!=i->second.end();j++)
	{
	  RangeAlign r1(0,j->getStart(),j->getEnd());
	  RangeAlign r2(0,i->first+j->getStart(),
			i->first+j->getEnd());

	  RangePair rp(r1,r2);
	  rp.setScore(r1.getLength());
	  sort_frag.insert(std::pair<unsigned long,RangePair>(rp.getScore(),rp));
	}
    }
}
//-------------------------------------------------------------------------
void HashAlign::print(unsigned min_size, std::ostream& out)
{
  out<<"Sequences:"<<std::endl;
  out<<sequence1.getDE()<<std::endl;
  out<<sequence2.getDE()<<std::endl;
  out<<"Alignment found:"<<std::endl;
  for(std::list<RangePair>::iterator i=rp_list.begin();
      i!=rp_list.end();i++)
    {
      if(i->getRangeQ().getLength()>=min_size
	 && i->getRangeS().getLength()>min_size)
	{
	  out<<i->getRangeQ().getStart()<<".."<<i->getRangeQ().getEnd()<<"\t"
	     <<i->getRangeS().getStart()<<".."<<i->getRangeS().getEnd()<<"\t"
	     <<i->getScore()<<std::endl;
	}
    } 
}
//-------------------------------------------------------------------------
void HashAlign::write(std::ostream& out)
{
  
  out<<"Sequences:"<<std::endl;
  out<<sequence1.getDE()<<std::endl;
  out<<sequence2.getDE()<<std::endl;
  int count=0;
  for(std::list<RangePair>::iterator i=rp_list.begin();
      i!=rp_list.end();i++)
    {
      out<<++count<<"\t"<<sequence1.getDE()<<"\t"
	 <<i->getRangeQ().getStart()<<"\t"<<i->getRangeQ().getEnd()<<std::endl;
      out<<count<<"\t"<<sequence2.getDE()<<"\t"
	 <<i->getRangeS().getStart()<<"\t"<<i->getRangeS().getEnd()<<std::endl;
    } 
}
//-------------------------------------------------------------------------
void HashAlign::fragAlign(double match,double mism, double gopen,
			   double gext, unsigned over,unsigned nfrag)
{
  FragAlign fragAligner(mism,gopen,gext,over);
  unsigned count=0;
  rp_list.clear();

  for(std::map<unsigned long,RangePair>::reverse_iterator i=sort_frag.rbegin();
      i!=sort_frag.rend();i++)
    {
      count++;

      rp_list.push_back(i->second);
      if(count>nfrag) break;
    }

  path=fragAligner.join(rp_list);
}
//-------------------------------------------------------------------------
// void HashAlign::plot(void)
// {
//   if(!h) return;
//   std::ostringstream cmd;
//   if(sequence1.length()!=0)
//     {
//       cmd<<"set xlabel \" "<<sequence1.getDE()<<" \" "<<std::endl
// 	 <<"set ylabel \" "<<sequence2.getDE()<<" \""<<std::endl
// 	 <<"set xrange [0:"<<sequence1.length()<<"]"<<std::endl
// 	 <<"set yrange [0:"<<sequence2.length()<<"]"<<std::endl;
//     }
//   else
//     {
//       unsigned minx=UINT_MAX,maxx=0,miny=UINT_MAX,maxy=0;
//       for(std::list<RangePairSet>::iterator i=path.begin();
// 	  i!=path.end();i++)
// 	for(std::list<RangePair>::iterator j=i->begin();
// 	    j!=i->end();j++)
// 	  {
// 	    minx=minx<std::min(j->getRangeQ().getStart(),j->getRangeQ().getEnd())?
// 	      minx:std::min(j->getRangeQ().getStart(),j->getRangeQ().getEnd());
// 	    miny=miny<std::min(j->getRangeS().getStart(),j->getRangeS().getEnd())?
// 	      miny:std::min(j->getRangeS().getStart(),j->getRangeS().getEnd());
// 	    maxx=maxx>std::max(j->getRangeQ().getStart(),j->getRangeQ().getEnd())?
// 	      maxx:std::max(j->getRangeQ().getStart(),j->getRangeQ().getEnd());
// 	    maxy=maxy>std::max(j->getRangeS().getStart(),j->getRangeS().getEnd())?
// 	      maxy:std::max(j->getRangeS().getStart(),j->getRangeS().getEnd());
// 	  }
//       cmd<<"set xrange ["<<minx<<":"<<maxx<<"]"<<std::endl
// 	 <<"set yrange ["<<miny<<":"<<maxy<<"]"<<std::endl;
//     }

//   std::ofstream fout1(datafile1);
	
//   for(std::list<RangePairSet>::iterator i=path.begin();
//       i!=path.end();i++)
//     for(std::list<RangePair>::iterator j=i->begin();
// 	j!=i->end();j++)
//       {
// 	fout1<<j->getRangeQ().getStart()
// 	    <<"\t"<<j->getRangeS().getStart()<<std::endl;
// 	fout1<<j->getRangeQ().getEnd()
// 	    <<"\t"<<j->getRangeS().getEnd()<<std::endl;
// 	fout1<<std::endl;
//       }
//   fout1.close();

//   std::ofstream fout2(datafile2);
	
//   for(std::list<RangePairSet>::iterator i=path.begin();
//       i!=path.end();i++)
//     for(std::list<RangePair>::iterator j=i->begin();
// 	j!=i->end();j++)
//       {
// 	std::list<RangePair>::iterator k=j;
// 	if(++k==i->end()) break;
// 	fout2<<j->getRangeQ().getEnd()
// 	    <<"\t"<<j->getRangeS().getEnd()<<std::endl;
// 	fout2<<k->getRangeQ().getStart()
// 	    <<"\t"<<k->getRangeS().getStart()<<std::endl;
// 	fout2<<std::endl;
//       }
//   fout2.close();

//   if(path.size())
//     {
//       cmd<<"plot \""<<datafile1<<"\" notitle with lines 1 1";
//       cmd<<", \""<<datafile2<<"\" notitle with lines 2 2";
//     }
//   gnuplot_cmd(h,(char*)cmd.str().c_str()) ;

// }






