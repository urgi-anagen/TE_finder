/**
 * \file BlastMatch.cpp
 */

#include <math.h>
#include <BlastMatch.h>

unsigned long BlastMatch::maxid=1;

void BlastMatch::reset(void)
{
  query_num=0;              
  subject_num=0;            

  query_strand='0';         
  subject_strand='0';       
  score=0;                  
  e_value=0;       
  identity=0;
  length=0;
  id=0;                     
  query_strt=subject_strt=0; 
  query_end=subject_end=0;  
}

std::istream& operator>>(std::istream& in, BlastMatch& align)
{
 
  char buff[256];
  align.reset();

  in.getline(buff,255);
  align.id=BlastMatch::maxid++;

  in.getline(buff,255);
  align.query_num=atol(buff);

  in.getline(buff,255);
  align.subject_num=atol(buff);

  in.getline(buff,255);
  align.query_strand=buff[0];

  in.getline(buff,255);
  align.subject_strand=buff[0];

  in.getline(buff,255);
  align.score=atol(buff);

  in.getline(buff,255);
  align.e_value=atof(buff);

  in.getline(buff,255);
  align.identity=atof(buff);

  in.getline(buff,255);
  align.length=atoi(buff);

  in.getline(buff,255);
  align.query_strt=atol(buff);

  in.getline(buff,255);
  align.query_end=atol(buff);

  in.getline(buff,255);
  align.subject_strt=atol(buff);

  in.getline(buff,255);
  align.subject_end=atol(buff);

  in.peek();

  return in;
}

std::istringstream& operator>>(std::istringstream& in, BlastMatch& align)
{
 
  char buff[256];
  align.reset();

  in.getline(buff,255);
  align.id=BlastMatch::maxid++;

  in.getline(buff,255);
  align.query_num=atol(buff);

  in.getline(buff,255);
  align.subject_num=atol(buff);

  in.getline(buff,255);
  align.query_strand=buff[0];

  in.getline(buff,255);
  align.subject_strand=buff[0];

  in.getline(buff,255);
  align.score=atol(buff);

  in.getline(buff,255);
  align.e_value=atof(buff);

  in.getline(buff,255);
  align.identity=atof(buff);

  in.getline(buff,255);
  align.length=atoi(buff);

  in.getline(buff,255);
  align.query_strt=atol(buff);

  in.getline(buff,255);
  align.query_end=atol(buff);

  in.getline(buff,255);
  align.subject_strt=atol(buff);

  in.getline(buff,255);
  align.subject_end=atol(buff);

  in.peek();

  return in;
}

std::ostream& operator<<(std::ostream& out, BlastMatch& align)
{
  if(align.query_num==0) return out;
  out<<align.id<<std::endl;
  out<<align.query_num<<std::endl;
  out<<align.subject_num<<std::endl;
  out<<align.query_strand<<std::endl;
  out<<align.subject_strand<<std::endl;
  out<<align.score<<std::endl;
  out<<align.e_value<<std::endl;
  out<<align.identity<<std::endl;
  out<<align.length<<std::endl;
  out<<align.query_strt<<std::endl;
  out<<align.query_end<<std::endl;
  out<<align.subject_strt<<std::endl;
  out<<align.subject_end<<std::endl;

  return out;
}

void BlastMatch::view(void) const
{
  if(query_num==0) return;
  std::cout<<id<<"\t";
  std::cout<<query_num<<"\t";
  std::cout<<subject_num<<"\t";
  std::cout<<query_strand<<"\t";
  std::cout<<subject_strand<<"\t";
  std::cout<<score<<"\t";
  std::cout<<e_value<<"\t";
  std::cout<<identity<<"\t";
  std::cout<<length<<"\t";
  std::cout<<query_strt<<"\t";
  std::cout<<query_end<<"\t";
  std::cout<<subject_strt<<"\t";
  std::cout<<subject_end<<std::endl;
}

void BlastMatch::readncbiblastfield(std::istream& in)
{

  char buff[256];
  reset();

  //query name
  in.getline(buff,255,'\t');   

  // le 1er champ de sortie sert a verifie si la sortie est vide ou non.
  if(in.eof())
    throw BlastMatch::Empty_Parser_output();

  query_num=atol(buff);

  //subject name
  in.getline(buff,255,'\t');
  subject_num=atol(buff);

  //identity %
  in.getline(buff,255,'\t');
  identity=atof(buff);

  //alignement length
  in.getline(buff,255,'\t');
  length=atoi(buff);

  //mismatches
  in.getline(buff,255,'\t');

  //gaps
  in.getline(buff,255,'\t');

  //query start
  in.getline(buff,255,'\t');
  query_strt=atol(buff);

  //query end
  in.getline(buff,255,'\t');
  query_end=atol(buff);

  //subject start
  in.getline(buff,255,'\t');
  subject_strt=atol(buff);
   
  //subject end
  in.getline(buff,255,'\t');
  subject_end=atol(buff);

  //e_value
  in.getline(buff,255,'\t');
  e_value=atof(buff);

  //score
  in.getline(buff,255);
  score=(unsigned long)rint(atof(buff));

  if(query_strt<query_end)
    {
      query_strand='0';
    }
  else
    {
      query_strand='1';
      unsigned long tmp=query_strt;
      query_strt=query_end;
      query_end=tmp;
    }

  if(subject_strt<subject_end)
    {
      subject_strand='0';
    }
  else
    {
      subject_strand='1';
      unsigned long tmp=subject_strt;
      subject_strt=subject_end;
      subject_end=tmp;
    }
  
  id=maxid++;

  in.peek();
}

void BlastMatch::readwublastfield(std::istream& in)
{

  char buff[256];
  reset();

  //query name
  in.getline(buff,255,'\t');   

  // le 1er champ de sortie sert a verifie si la sortie est vide ou non.
  if(in.eof())
    throw BlastMatch::Empty_Parser_output();

  query_num=atol(buff);

  //subject name
  in.getline(buff,255,'\t');
  subject_num=atol(buff);

  //e_value
  in.getline(buff,255,'\t');
  e_value=atof(buff);


  //number of score for e-value
  in.getline(buff,255,'\t');

  //score
  in.getline(buff,255,'\t');
  score=(unsigned long)rint(atof(buff));

  //raw score
  in.getline(buff,255,'\t');

  //alignement length
  in.getline(buff,255,'\t');
  length=atoi(buff);

  //matches
  in.getline(buff,255,'\t');
  unsigned nmatch=atoi(buff);

  //positives matches
  in.getline(buff,255,'\t');

  //mismatches (negative matches)
  in.getline(buff,255,'\t');

  //identity %
  in.getline(buff,255,'\t');

  //sim %
  in.getline(buff,255,'\t');

  //number of gaps in query
  in.getline(buff,255,'\t');
  //unsigned ngapq=atoi(buff);

  //total gaps length in query
  in.getline(buff,255,'\t');
  unsigned lgapq=atoi(buff);

  //number of gaps in subject
  in.getline(buff,255,'\t');
  //unsigned ngaps=atoi(buff);

  //total gaps length in subject
  in.getline(buff,255,'\t');
  unsigned lgaps=atoi(buff);

  identity=((float)nmatch/(length-lgapq-lgaps))*100;

  //frame in query
  in.getline(buff,255,'\t');

  //query start
  in.getline(buff,255,'\t');
  query_strt=atol(buff);

  //query end
  in.getline(buff,255,'\t');
  query_end=atol(buff);

  //frame in subject
  in.getline(buff,255,'\t');

  //subject start
  in.getline(buff,255,'\t');
  subject_strt=atol(buff);
   
  //subject end
  in.getline(buff,255);
  subject_end=atol(buff);

  if(query_strt<query_end)
    {
      query_strand='0';
    }
  else
    {
      query_strand='1';
      unsigned long tmp=query_strt;
      query_strt=query_end;
      query_end=tmp;
    }

  if(subject_strt<subject_end)
    {
      subject_strand='0';
    }
  else
    {
      subject_strand='1';
      unsigned long tmp=subject_strt;
      subject_strt=subject_end;
      subject_end=tmp;
    }
  
  id=maxid++;

  in.peek();
}

void BlastMatch::readlst(std::istream& in)
{
  in.read((char*)&id, sizeof(id));
  in.read((char*)&query_num, sizeof(query_num));
  in.read((char*)&subject_num, sizeof(subject_num));
  in.read((char*)&query_strand, sizeof(query_strand));
  in.read((char*)&subject_strand, sizeof(subject_strand));
  in.read((char*)&score, sizeof(score));
  in.read((char*)&e_value, sizeof(e_value));
  in.read((char*)&identity, sizeof(identity));
  in.read((char*)&length, sizeof(length));
  in.read((char*)&query_strt, sizeof(query_strt));
  in.read((char*)&query_end, sizeof(query_end));
  in.read((char*)&subject_strt, sizeof(subject_strt));
  in.read((char*)&subject_end, sizeof(subject_end));
  in.peek();
}

void BlastMatch::readoldlst(std::istream& in)
{
  in.read((char*)&id, sizeof(id));
  in.read((char*)&query_num, sizeof(query_num));
  in.read((char*)&subject_num, sizeof(subject_num));
  in.read((char*)&query_strand, sizeof(query_strand));
  in.read((char*)&subject_strand, sizeof(subject_strand));
  in.read((char*)&score, sizeof(score));
  in.read((char*)&e_value, sizeof(e_value));
  in.read((char*)&identity, sizeof(identity));
  in.read((char*)&query_strt, sizeof(query_strt));
  in.read((char*)&query_end, sizeof(query_end));
  in.read((char*)&subject_strt, sizeof(subject_strt));
  in.read((char*)&subject_end, sizeof(subject_end));
  in.peek();
}

void BlastMatch::writelst(std::ostream& out)
{
  if(query_num==0) return;
  out.write((char*)&id, sizeof(id));
  out.write((char*)&query_num, sizeof(query_num));
  out.write((char*)&subject_num, sizeof(subject_num));
  out.write((char*)&query_strand, sizeof(query_strand));
  out.write((char*)&subject_strand, sizeof(subject_strand));
  out.write((char*)&score, sizeof(score));
  out.write((char*)&e_value, sizeof(e_value));
  out.write((char*)&identity, sizeof(identity));
  out.write((char*)&length, sizeof(length));
  out.write((char*)&query_strt, sizeof(query_strt));
  out.write((char*)&query_end, sizeof(query_end));
  out.write((char*)&subject_strt, sizeof(subject_strt));
  out.write((char*)&subject_end, sizeof(subject_end)); 
}
