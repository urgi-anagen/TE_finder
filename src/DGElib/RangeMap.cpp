/*
 *
 * RangeMap.cpp
 *
 */

#include <stdio.h>
#include <math.h>
#include "RangeMap.h"

//-------------------------------------------------------------------
void RangeMap::load(const SDGString& filename)
{
  std::ifstream fin((const char*)filename.start());

  if( ! fin.is_open() )
    {
      std::cerr<<"Error opening file "<<filename<<std::endl;
      exit(1);
    }
  while(!fin.eof())
    {
      RangeSeq r;
      fin>>r;
      if(!r.empty()) add(r);
      fin.peek();
    }
  fin.close();
  sort();
}
//-------------------------------------------------------------------
void RangeMap::save(const SDGString& filename)
{
  std::ofstream fout((const char*)filename.start());

  for(RangeMap::iterator c=begin();
	c!=end();c++)
      {
	for(std::list<RangeSeq>::iterator i=c->second.begin();
	    i!=c->second.end();i++)
	  {
		fout<<(*i)<<std::endl;
	  }
      }
  fout.close();
}
//-------------------------------------------------------------------
void RangeMap::saveSet(const SDGString& filename)
{
  std::ofstream fout((const char*)filename.start());

  unsigned id=0;
  for(RangeMap::iterator c=begin();
	c!=end();c++)
      {
	for(std::list<RangeSeq>::iterator i=c->second.begin();
	    i!=c->second.end();i++)
	  {
	    i->write_rangeSet(fout,++id);
	  }
      }
  fout.close();
}
//-------------------------------------------------------------------
void RangeMap::sort(void)
{
  for(RangeMap::iterator c=begin();
	c!=end();c++)
	  c->second.sort();
}
//-------------------------------------------------------------------
void RangeMap::view(void)
{
  int count=0;
  for(RangeMap::iterator c=begin(); c!=end();c++)
    {
      ulong prev=0;
      for(std::list<RangeSeq>::iterator i=c->second.begin();
	  i!=c->second.end();i++)
	{
	  long is=(*i).getStart();
	  long ie=(*i).getEnd();
	  long long dist=(is<ie?is:ie)-prev;
	  std::cout<<"#"<<++count<<"\t"
	      <<(*i)
	      <<"\tlen="<<(unsigned)std::abs(ie-is+1)<<"\t"
	      <<"dist="<<dist;
	  if(is<ie)
	    std::cout<<"\t+";
	  else
	    std::cout<<"\t-";
	  if(prev==0) std::cout<<"\t(first of sequence)";
	  std::cout<<std::endl;
	  prev=(is>ie?is:ie);
	}
    }
}
//-------------------------------------------------------------------
unsigned RangeMap::size(void)
{
  unsigned size=0;
  for(RangeMap::iterator c=begin(); c!=end();c++)
    {
      for(std::list<RangeSeq>::iterator i=c->second.begin();
	  i!=c->second.end();i++)
		{
		  long is=(*i).getStart();
		  long ie=(*i).getEnd();
		  size+=(unsigned)std::abs(ie-is+1);
		}
    }
  return size;
}
//-------------------------------------------------------------------
void RangeMap::cut(unsigned int length, unsigned int over )
{
  for(RangeMap::iterator c=begin(); c!=end();c++)
    {
      for(std::list<RangeSeq>::iterator i=c->second.begin();
	  i!=c->second.end();i++)
	{
	  long is=(*i).getStart();
	  long ie=(*i).getEnd();
	  if((unsigned)std::abs(ie-is+1)>length)
	    {
	      (*i).setStart(is+length-over);
	      RangeSeq new_range((*i).getName(),(*i).getChr(),
				    is,is+length-1);
	      c->second.insert(i,new_range);
	      i--;
	    }
	  else
	    {
	      (*i).setName((*i).getName());
	    }
	}
    }
}
//-------------------------------------------------------------------
void RangeMap::selectSrcSeq(const SDGString& outfname,const SDGBioSeqDB& db)
{
  SDGFastaOstream out(outfname);

  for(SDGBioSeqDB::const_iterator db_it=db.begin();
      db_it!=db.end();db_it++)
    {
      SDGBioSeq chr=(*db_it);
      std::string h_seq=(const char*)(chr.getDE().start());
      RangeMap::iterator c=find(h_seq);
      if(c!=end())
	out<<chr;
    }
}
//-------------------------------------------------------------------
void RangeMap::extend(unsigned size_flank)
{
  for(RangeMap::iterator c=begin(); c!=end();c++)
    {
      for(std::list<RangeSeq>::iterator i=c->second.begin();
	  i!=c->second.end();i++)
	    {
	      if(i->isPlusStrand())
		{
		  i->setStart(size_flank<(*i).getStart()?(*i).getStart()-size_flank:1);
		  i->setEnd(i->getEnd()+size_flank);
		}
	      else
		{
		  i->setStart(i->getStart()+size_flank);
		  i->setEnd(size_flank<(*i).getEnd()?(*i).getEnd()-size_flank:1);
		}
	    }
    }
}
//-------------------------------------------------------------------
void RangeMap::writeSeq(const SDGString& outfname,const SDGBioSeqDB& db) const
{
  SDGFastaOstream out(outfname);
  int count_skip=0;

  for(SDGBioSeqDB::const_iterator db_it=db.begin();
      db_it!=db.end();db_it++)
    {
      SDGBioSeq chr=(*db_it);
      unsigned chr_len=chr.length();
      std::string h_seq=(const char*)(chr.getDE().start());
      RangeMap::const_iterator c=find(h_seq);
      if(c!=end())
	{
	  for(std::list<RangeSeq>::const_iterator i=c->second.begin();
	      i!=c->second.end();i++)
	    {
	      if(i->getLength()<14)// temporaire !!??
		{
		  count_skip++;
		  continue;
		}
	      long is=i->getStart();
	      long ie=i->getEnd();
	      if(is<=ie)
		{
		  if(is>=chr_len)
		    continue;
		}
	      else
		{
		  if(ie>=chr_len)
		    continue;
		}
	      SDGBioSeq s=newSDGMemBioSeq("");
	      std::ostringstream name;
	      if(i->getRangeSet().size()<=1)
		{
		  if(is<=ie)
		    {
		      s=chr.subseq(is-1,ie-is+1);
		    }
		  else
		    {
		      s=chr.subseq(ie-1,is-ie+1);
		      s=s.complement();
		    }
		  name<<i->getName()<<" "<<i->getChr()
		      <<" {Fragment} "<<is<<".."<<ie;
		}
	      else
		{
		  SDGBioSeq sr=newSDGMemBioSeq("");
		  name<<i->getName()<<" "<<i->getChr()
		      <<" {Fragment} ";
		  bool first=true;
		  for(std::list<Range>::const_iterator r= i->getRangeSet().begin();
		      r!=i->getRangeSet().end();r++)
		    {
		      ulong rs=r->getStart();
		      ulong re=r->getEnd();
		      if(r->isPlusStrand())
			sr=chr.subseq(rs-1,re-rs+1);
		      else
			{
			  sr=chr.subseq(re-1,rs-re+1);
			  sr=sr.complement();
			}
		      s+=sr;
		      if(!first) name<<",";
		      else first=false;
		      name<<rs<<".."<<re;
		    }
		}
	      s.setDE(name.str());
	      out<<s;
	    }
	}
    }
  if(count_skip>0)
	  std::cout<<count_skip<<" skipped!"<<std::endl;
}
//-------------------------------------------------------------------
void RangeMap::writeCutSeq( const SDGString& outfname, const SDGBioSeqDB& db, int verbose )
{
  SDGFastaOstream out(outfname);
  int countTreated=0;
  int count_skip=0;

  for(SDGBioSeqDB::const_iterator db_it=db.begin();
      db_it!=db.end();db_it++)
    {
      SDGBioSeq chr=(*db_it);
      std::string h_seq=(const char*)(chr.getDE().start());
      RangeMap::const_iterator c=find(h_seq);
      if(c!=end())
	{
	  for(std::list<RangeSeq>::const_iterator i=c->second.begin();
	      i!=c->second.end();i++)
	    {
            long is = i->getStart();
            long ie = i->getEnd();
            if (i->getLength() < 14)// temporaire !!
            {
                count_skip++;
                if(verbose>0)
                    std::cout << "Squence length < 14 .. skip sequence:" << i->getName() << " " << i->getChr()
                              << " {Cut} " << is << ".." << ie << std::endl;
                continue;
            }

	      countTreated ++;
	      if(verbose>0)
	      {
	    	  std::cout<<"   working on range "<<countTreated<<"/"<<countRange
	    	  <<" length="<<(long)std::abs(is-ie)+1<<std::endl;
	    	  std::cout<<*i<<std::endl<<std::flush;
	      }
	      SDGBioSeq s=newSDGMemBioSeq("");
	      if(is<=ie)
		{
		  s=chr.subseq(is-1,ie-is+1);
		}
	      else
		{
		  s=chr.subseq(ie-1,is-ie+1);
		  s=s.complement();
		}
	      std::ostringstream name;
	      name<<countTreated<<" "<<i->getName()<<" "<<i->getChr()
		  <<" {Cut} "<<is<<".."<<ie;
	      s.setDE(name.str());
	      out<<s;
	    }
	}
    }
  if(count_skip>0)
	  if(verbose>0)
		  std::cout<<count_skip<<" skipped!"<<std::endl;
}
//-------------------------------------------------------------------
void RangeMap::writeFlank53Seq(const SDGString& outfname,const SDGBioSeqDB& db,unsigned len)
{
  SDGFastaOstream out(outfname);
  int countTreated=0;
  int count_skip=0;

  for(SDGBioSeqDB::const_iterator db_it=db.begin();
      db_it!=db.end();db_it++)
    {
      SDGBioSeq chr=(*db_it);
      std::string h_seq=(const char*)(chr.getDE().start());
      RangeMap::const_iterator c=find(h_seq);
      if(c!=end())
	{
	  for(std::list<RangeSeq>::const_iterator i=c->second.begin();
	      i!=c->second.end();i++)
	    {
	      if(i->getLength()<14)// temporaire !!
		{
		  count_skip++;
		  continue;
		}

	      std::cout<<"   working on range "<<++countTreated<<"/"<<countRange
		  <<" length="<<i->getLength()<<std::endl;
	      std::cout<<*i<<std::endl<<std::flush;

	      SDGBioSeq s1,s2,s3;
	      ulong is1,ie1,is2,ie2;
	      if(i->getStart()<i->getEnd())
		{
		  // 5' flank

		  is1=len<i->getStart()?i->getStart()-len:1;
		  ie1=i->getStart();
		  s1=chr.subseq(is1-1,ie1-is1+1);

		  // 3' flank
		  is2=i->getEnd();
		  ie2=i->getEnd()+len;
		  s2=chr.subseq(is2-1,ie2-is2+1);

		  s3=newSDGMemBioSeq((SDGString)(s1.toString()+s2.toString()));
		}
	      else
		{
		  // 3' flank
		  is1=len<i->getEnd()?i->getEnd()-len:1;
		  ie1=i->getEnd();
		  s1=chr.subseq(is1-1,ie1-is1+1);

		  // 5' flank
		  is2=i->getStart();
		  ie2=i->getStart()+len;
		  s2=chr.subseq(is2-1,ie2-is2+1);

		  s3=newSDGMemBioSeq(s1.toString()+s2.toString());
		  s3=s3.complement();
		}


	      //write seq
	      char start1[256],end1[256];
	      sprintf(start1,"%ld",is1);
	      sprintf(end1,"%ld",ie1);
	      char start2[256],end2[256];
	      sprintf(start2,"%ld",is2);
	      sprintf(end2,"%ld",ie2);
	      std::string name=i->getName()+" "+i->getChr()+
		" "+start1+".."+end1+","+start2+".."+end2;
	      s3.setDE(SDGString(name.c_str()));
	      out<<s3;
	    }
	}
    }
  if(count_skip>0)
	  std::cout<<count_skip<<" skipped!"<<std::endl;
}
//-------------------------------------------------------------------
void RangeMap::writeFlank5Seq(const SDGString& outfname,const SDGBioSeqDB& db,unsigned len)
{
  SDGFastaOstream out(outfname);
  int countTreated=0;
  int count_skip=0;

  for(SDGBioSeqDB::const_iterator db_it=db.begin();
      db_it!=db.end();db_it++)
    {
      SDGBioSeq chr=(*db_it);
      std::string h_seq=(const char*)(chr.getDE().start());
      RangeMap::const_iterator c=find(h_seq);
      if(c!=end())
	{
	  for(std::list<RangeSeq>::const_iterator i=c->second.begin();
	      i!=c->second.end();i++)
	    {
	      if(i->getLength()<14)// temporaire !!
		{
		  count_skip++;
		  continue;
		}

	      std::cout<<"   working on range "<<++countTreated<<"/"<<countRange
		  <<" length="<<i->getLength()<<std::endl;
	      std::cout<<*i<<std::endl<<std::flush;

	      SDGBioSeq s;
	      ulong is,ie;
	      if(i->getStart()<i->getEnd())
		{
		  // 5' flank
		  is=len<i->getStart()?i->getStart()-len:1;
		  ie=i->getStart();
		  s=chr.subseq(is-1,ie-is+1);
		}
	      else
		{
		  // 5' flank
		  is=i->getStart();
		  ie=i->getStart()+len;
		  s=chr.subseq(is-1,ie-is+1);
		  s=s.complement();
		}


	      //write seq
	      char start[256],end[256];
	      sprintf(start,"%ld",is);
	      sprintf(end,"%ld",ie);
	      std::string name=i->getName()+" "+i->getChr()+
		" 5' flank "+start+".."+end;
	      s.setDE(SDGString(name.c_str()));
	      out<<s;
	    }
	}
    }
  if(count_skip>0)
	  std::cout<<count_skip<<" skipped!"<<std::endl;
}
//-------------------------------------------------------------------
void RangeMap::writeFlank3Seq(const SDGString& outfname,const SDGBioSeqDB& db,unsigned len)
{
  SDGFastaOstream out(outfname);
  int countTreated=0;
  int count_skip=0;

  for(SDGBioSeqDB::const_iterator db_it=db.begin();
      db_it!=db.end();db_it++)
    {
      SDGBioSeq chr=(*db_it);
      std::string h_seq=(const char*)(chr.getDE().start());
      RangeMap::const_iterator c=find(h_seq);
      if(c!=end())
	{
	  for(std::list<RangeSeq>::const_iterator i=c->second.begin();
	      i!=c->second.end();i++)
	    {
	      if(i->getLength()<14)// temporaire !!
		{
		  count_skip++;
		  continue;
		}

	      std::cout<<"   working on range "<<++countTreated<<"/"<<countRange
		  <<" length="<<i->getLength()<<std::endl;
	      std::cout<<*i<<std::endl<<std::flush;

	      SDGBioSeq s;
	      ulong is,ie;
	      if(i->getStart()<i->getEnd())
		{
		  // 3' flank
		  is=i->getEnd();
		  ie=i->getEnd()+len;
		  s=chr.subseq(is-1,ie-is+1);
		}
	      else
		{
		  // 3' flank
		  is=len<i->getEnd()?i->getEnd()-len:1;
		  ie=i->getEnd();
		  s=chr.subseq(is-1,ie-is+1);
		  s=s.complement();
		}


	      //write seq
	      char start[256],end[256];
	      sprintf(start,"%ld",is);
	      sprintf(end,"%ld",ie);
	      std::string name=i->getName()+" "+i->getChr()+
		" 3' flank "+start+".."+end;
	      s.setDE(SDGString(name.c_str()));
	      out<<s;
	    }
	}
    }
  if(count_skip>0)
	  std::cout<<count_skip<<" skipped!"<<std::endl;
}
//-------------------------------------------------------------------
void RangeMap::merge(void)
{
  bool isReajusted=true;
  while(isReajusted==true)
    {
      isReajusted=false;
      for(RangeMap::iterator c=begin(); c!=end();c++)
	{
	  for(std::list<RangeSeq>::iterator i=c->second.begin();
	      i!=c->second.end();i++)
	    {
	      std::list<RangeSeq>::iterator j=i;
	      j++;
	      while(j!=c->second.end() && i->overlap(*j))
		{
		  i->merge(*j);
		  isReajusted=true;
		  j=c->second.erase(j);
		  countRange--;
		}
	    }
	}
    }
}
//-------------------------------------------------------------------
void RangeMap::diff( RangeMap& m, int verbose ) {
    for (RangeMap::iterator c = begin(); c != end(); c++) {
        for (std::list<RangeSeq>::iterator i = c->second.begin();
             i != c->second.end(); i++) {
            std::list<RangeSeq>::iterator j;
            for (j = m[c->first].begin(); i != c->second.end()
                                          && j != m[c->first].end(); j++)
                if (i->overlap(*j)) {
                    if (verbose > 0) {
                        std::cout << "\nsubject:\t" << (*i) << std::endl;
                        std::cout << "query:\t" << (*j) << std::endl;
                    }
                    RangeSeq r = i->diff(*j);
                    if (i->empty()) {
                        i = c->second.erase(i);
                        i--;
                        if (verbose > 0)
                            std::cout << "==>erase subject!" << std::endl;
                        break;
                    }
                    if (verbose > 0)
                        std::cout << "==>change subject:" << (*i) << std::endl;
                    if (!r.empty()) {
                        if (verbose > 0)
                            std::cout << "==>new subject:" << r << std::endl;
                        c->second.push_back(r);
                    }
                }
        }
        c->second.sort();
    } // for
}
//-------------------------------------------------------------------
void RangeMap::selectInclude(RangeMap& m, RangeMap& mapout)
{
  for(RangeMap::iterator c=begin(); c!=end();c++)
    {
      for(std::list<RangeSeq>::iterator i=c->second.begin();
	  i!=c->second.end();i++)
	{
	  std::list<RangeSeq>::iterator j;
	  for( j=m[c->first].begin();  i!=c->second.end()
		 && j!=m[c->first].end();j++)
	      if(i->isIncluded(*j))
		{
		  mapout.add(*i);
		  break;
		}
	}
    }
}
//-------------------------------------------------------------------
void RangeMap::selectOverlap(RangeMap& m, RangeMap& mapout)
{
  for(RangeMap::iterator c=begin(); c!=end();c++)
    {
      for(std::list<RangeSeq>::iterator i=c->second.begin();
	  i!=c->second.end();i++)
	{
	  std::list<RangeSeq>::iterator j;
	  for( j=m[c->first].begin();  i!=c->second.end()
		 && j!=m[c->first].end();j++)
	    if(i->overlap(*j))
	      {
		mapout.add(*i);
		break;
	      }
	}
    }
}
//-------------------------------------------------------------------
void RangeMap::entropfilt(const SDGString& db, double thres)
{
  SDGFastaIstream in(db);
  int countTreated=0;
  int count_skip=0;
  int countErase=0;

  SDGBioSeq chr;
  while(in)
    {
      in>>chr;
      std::string h_seq=(const char*)(chr.getDE().start());
      RangeMap::iterator c=find(h_seq);
      if(c!=end())
	{
	  for(std::list<RangeSeq>::iterator
		i=c->second.begin();
	      i!=c->second.end();i++)
	    {
	      if(i->getLength()==1)// temporaire !!
		{
		  count_skip++;
		  continue;
		}
	      long is=i->getStart();
	      long ie=i->getEnd();
	      std::cout<<"   working on range "<<++countTreated<<"/"<<countRange
		  <<" length="<<(long)std::abs(is-ie)+1<<std::endl;
	      std::cout<<*i<<std::endl<<std::flush;
	      SDGBioSeq s=newSDGMemBioSeq("");
	      if(is<=ie)
		{
		  s=chr.subseq(is-1,ie-is+1);
		}
	      else
		{
		  s=chr.subseq(ie-1,is-ie+1);
		  s=s.complement();
		}
	      unsigned len=s.length();
	      int countA=0;
	      int countT=0;
	      int countG=0;
	      int countC=0;
	      for(unsigned p=0;p<len;p++)
		{
		  switch(s.charAt(p))
		    {
		    case 'A' : countA++; break;
		    case 'T' : countT++; break;
		    case 'G' : countG++; break;
		    case 'C' : countC++; break;
		    }
		}
	      countA++;
	      countT++;
	      countG++;
	      countC++;
	      double sum=countA+countT+countG+countC;
	      double pA=countA/sum;
	      double pT=countT/sum;
	      double pG=countG/sum;
	      double pC=countC/sum;
	      double ent=(-1)*(pA * log(pA)/log(4)
			   + pT * log(pT)/log(4)
			   + pG * log(pG)/log(4)
			   + pC * log(pC)/log(4));
	      std::cout<<"\tentropy="<<ent;
	      if(ent<thres)
		{
		  std::cout<<" low entropy!";
		  std::list<RangeSeq>::iterator tmp=i;
		  i--;
		  c->second.erase(tmp);
		  countErase++;
		}
	      std::cout<<std::endl;
	    }
	}
    }
  std::cout<<"Erased #:"<<countErase<<std::endl;
}
