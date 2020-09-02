/*
 * \file RangePair.cpp
 */

#include <vector>
#include "RangePair.h"
#include "SDGString.h"

const RangePair::Less RangePair::less;
const RangePair::Greater RangePair::greater;
const RangePair::GreaterScore RangePair::greaterScore;
const RangePair::GreaterScoreIdLenCoord RangePair::greaterScoreIdLenCoord;
const RangePair::LessIdentity RangePair::lessIdentity;
const RangePair::GreaterLengthQ RangePair::greaterLengthQ;
const RangePair::GreaterLengthIdent RangePair::greaterLengthIdent;
const RangePair::StrictLess RangePair::strictLess;


RangePair::RangePair( SDGString line )
{
	set( line );
}

void RangePair::set( SDGString line )
{
	id=0;
	std::vector<SDGString> tokens = SDGString::tokenize( line );
	if(tokens.size() != 9)
	{
		if(tokens.size() != 6)
		{
			std::cerr<<"Info on current RangePair created: "<<tokens[0]<<"\t"<<tokens[1]<<"\t"<<tokens[2]<<"\t"<<tokens[3]<<"\t"<<tokens[4]<<"\t"<<tokens[5]<<"\t"<<std::endl<<std::flush;
		}
		std::cerr<<"ERROR: not enough data in line to set RangePair"<<std::endl<<std::flush;
		std::exit( EXIT_FAILURE );
	}
	first.set( tokens[0], 0, atoi(tokens[1]), atoi(tokens[2]) );
	second.set( tokens[3], 0, atoi(tokens[4]), atoi(tokens[5]) );
	e_value = atof(tokens[6]);
	score = atoi(tokens[7]);
	identity = atof(tokens[8]);
	length = std::max(std::abs((double)(first.getEnd())-first.getStart())+1,std::abs((double)(second.getEnd())-second.getStart())+1);
}

bool operator==( const RangePair &rp1, const RangePair &rp2 )
{
	return( rp1.first == rp2.first
			&& rp1.second == rp2.second
			&& rp1.e_value == rp2.e_value
			&& rp1.score == rp2.score
			&& rp1.identity == rp2.identity
			&& rp1.length == rp2.length);
}

RangePair::RangePair(BlastMatch al)
{
  score=al.getScore();
  e_value=al.getE_value();
  id=al.getId();
  identity=al.getIdentity();
  length=al.getLength();

  if(al.getQuery_strand()=='0')
    first.set(al.getQuery_num(),
	       al.getQuery_strt(),al.getQuery_end());
  else
    first.set(al.getQuery_num(),
	       al.getQuery_end(),al.getQuery_strt());

  if(al.getSubject_strand()=='0')
    second.set(al.getSubject_num(),
	       al.getSubject_strt(),al.getSubject_end());
  else
    second.set(al.getSubject_num(),
	       al.getSubject_end(),al.getSubject_strt());
  setStrand();
};

BlastMatch RangePair::toAlign(void)
{
  BlastMatch al;

  al.setQuery_num(first.getNumChr());
  al.setSubject_num(second.getNumChr());

  al.setQuery_strand(first.isPlusStrand()?'0':'1');
  al.setSubject_strand(second.isPlusStrand()?'0':'1');

  al.setQuery_strt(std::min(first.getStart(),first.getEnd()));
  al.setSubject_strt(std::min(second.getStart(),second.getEnd()));

  al.setQuery_end(std::max(first.getStart(),first.getEnd()));
  al.setSubject_end(std::max(second.getStart(),second.getEnd()));

  al.setE_value(e_value);
  al.setScore(score);
  al.setId(id);
  al.setIdentity(identity);
  al.setLength(length);

  return al;
};

void RangePair::readlst(std::istream& in)
{
  in.read((char*)&id, sizeof(id));

  unsigned long query_num;
  in.read((char*)&query_num, sizeof(query_num));
  unsigned long subject_num;
  in.read((char*)&subject_num, sizeof(subject_num));

  char query_strand;
  in.read((char*)&query_strand, sizeof(query_strand));
  char subject_strand;
  in.read((char*)&subject_strand, sizeof(subject_strand));

  in.read((char*)&score, sizeof(score));
  in.read((char*)&e_value, sizeof(e_value));
  in.read((char*)&identity, sizeof(identity));
  in.read((char*)&length, sizeof(length));

  unsigned long query_strt, query_end;
  in.read((char*)&query_strt, sizeof(query_strt));
  in.read((char*)&query_end, sizeof(query_end));
  if(query_strand=='0')
    first.set(query_num,
	       query_strt,query_end);
  else
    first.set(query_num,
	       query_end,query_strt);

  unsigned long subject_strt,subject_end;
  in.read((char*)&subject_strt, sizeof(subject_strt));
  in.read((char*)&subject_end, sizeof(subject_end));
  if(subject_strand=='0')
    second.set(subject_num,
	       subject_strt,subject_end);
  else
    second.set(subject_num,
	       subject_end,subject_strt);

  setStrand();
  in.peek();
}

void RangePair::readReputer(std::istream& in)
{
  char buff[256];

  in.getline(buff,255);
  if(buff[0]=='#') return;

  std::istringstream line(buff);
  unsigned length1=0;
  line>>length1;

  unsigned start1=0;
  line>>start1;

  char type='?';
  line>>type;

  unsigned length2=0;
  line>>length2;

  unsigned start2=0;
  line>>start2;

  int mismh=0;
  line>>mismh;


  line>>e_value;

  first.setStart(start1+1);
  first.setEnd(start1+length1);

  if(type=='F')
    {
      second.setStart(start2+1);
      second.setEnd(start2+length2);
    }
  else
    if(type=='P')
    {
      second.setStart(start2+length2);
      second.setEnd(start2+1);
    }

  length=std::max(length1,length2);
  score=length-abs(mismh);
  identity=(double)(length-abs(mismh))/length;

  first.setNumChr(1);
  second.setNumChr(1);
  setStrand();
}

void RangePair::readtxt(std::istream& in)
{
  char buff[1024];

  SDGString qname,sname;
  long qstart,qend,sstart,send;

  in.getline(buff,1023,'\t');
  qname=buff;

  in.getline(buff,1023,'\t');
  qstart=atol(buff);

  in.getline(buff,1023,'\t');
  qend=atol(buff);

  in.getline(buff,1023,'\t');
  sname=buff;

  in.getline(buff,1023,'\t');
  sstart=atol(buff);

  in.getline(buff,1023,'\t');
  send=atol(buff);

  in.getline(buff,1023,'\t');
  e_value=atof(buff);

  in.getline(buff,1023,'\t');
  score=atol(buff);

  in.getline(buff,1023,'\n');
  identity=atof(buff);


  first.set(qname,0,qstart,qend);
  second.set(sname,0,sstart,send);

  length=std::max(std::abs(qend-qstart)+1,std::abs(send-sstart)+1);

  setStrand();
}

void RangePair::writetxt(std::ostream& out) const
{
  out<<first.getNameSeq()<<"\t";
  out<<first.getStart()<<"\t";
  out<<first.getEnd()<<"\t";
  out<<second.getNameSeq()<<"\t";
  out<<second.getStart()<<"\t";
  out<<second.getEnd()<<"\t";
  out<<e_value<<"\t";
  out<<score<<"\t";
  out<<identity<<"\n";
}

void RangePair::merge(RangePair& r)
{
  first.merge(r.first);
  second.merge(r.second);
  score+=r.score;
  identity=(((identity/100)*length+(r.identity/100)*r.length)/(length+r.length))*100;
  e_value=std::min(e_value,r.e_value);
  length=std::abs((int)(first.getStart()-first.getEnd()));
}

void RangePair::merge(const RangePair& r)
{
  first.merge(r.first);
  second.merge(r.second);
  score+=r.score;
  identity=(((identity/100)*length+(r.identity/100)*r.length)/(length+r.length))*100;
  e_value=std::min(e_value,r.e_value);
  length=std::abs((int)(first.getStart()-first.getEnd()));
}

RangePair RangePair::diffQ(const RangePair& r)
{
      unsigned qs=first.getStart();
      unsigned qe=first.getEnd();
      RangeAlign ra=first.diff(r.first);
      RangePair new_r;

      if (!ra.empty())
		{
		  new_r=*this;
		  new_r.getRangeQ()=ra;
		}

      reComputeSubjectCoords(new_r,qs,qe);

      if (!ra.empty())
		{
		  unsigned s=(unsigned)floor(score*new_r.computeSizeRatioOnQuery(qs,qe)+0.5);
		  s=s>0?s:0;
		  new_r.setScore((unsigned long)s);

		}
      double ratio=computeSizeRatioOnQuery(qs,qe);
      length=(unsigned)floor(length*(ratio)+0.5);
      long s=(long)floor(score*(ratio)+0.5);
      s=s>0?s:0;
      score=(long)s;
      return new_r;
};

void RangePair::reComputeSubjectCoords(RangePair& new_r, unsigned qs, unsigned qe)
{
   double ratio_self=computeSizeRatioOnQuery(qs,qe);
   double ratio_new=new_r.computeSizeRatioOnQuery(qs,qe);

   unsigned size_subj_self=(unsigned)floor(ratio_self*second.getLength()+0.5);
   unsigned size_rm_self=second.getLength()-size_subj_self;

   unsigned size_subj_new=(unsigned)floor(ratio_new*new_r.getRangeS().getLength()+0.5);
   unsigned size_rm_new=new_r.getRangeS().getLength()-size_subj_new;

   unsigned long localStart = 0;
   unsigned long localEnd = 0;

	if(isQueryStartHasMoved(qs))
	{
	  if(isPlusStrand())
	    {
	      localStart = second.getStart() + size_rm_self;
	      localEnd = new_r.getRangeS().getEnd()-size_rm_new;
	    }
	  else
	    {
	      localStart = second.getStart()-size_rm_self;
	      localEnd = new_r.getRangeS().getEnd()+size_rm_new;
	     
	    }
	  second.setStart(localStart);
	  new_r.getRangeS().setEnd(localEnd);
	}

	if(isQueryEndHasMoved(qe))
	{

	  if(isPlusStrand())
	    {
	      localEnd = second.getEnd()-size_rm_self;
	      localStart = new_r.getRangeS().getStart()+size_rm_new;
	    }
	  else
	    {
	      localEnd = second.getEnd()+size_rm_self;
	      localStart = new_r.getRangeS().getStart()-size_rm_new;
	    }
	  second.setEnd(localEnd);
	  new_r.getRangeS().setStart(localStart);
	}
};

