/*
 *
 * BlastMatch.h
 *
 */

#ifndef BLASTMATCH_H
#define BLASTMATCH_H

#include <stdlib.h>
#include <iostream>
#include <sstream>

class BlastMatch 
{
 private:
  unsigned long query_num;     // sequence number in query bank
  unsigned long subject_num;   // sequence number in subject bank
  char query_strand;           // Strand of query
  char subject_strand;         // Strand of subject 
  unsigned long score;         // Score of the alignment 
  double e_value;         
  double identity;         
  unsigned length;         
  unsigned long id;           // Number of the match
  unsigned long query_strt;   // start of query
  unsigned long query_end;    // end of query 
  unsigned long subject_strt; // start of subject
  unsigned long subject_end;  // end of subject

  
 public:
  class Empty_Parser_output
    {
    public:
       Empty_Parser_output (void){};
    };
  
  static unsigned long maxid;     // variable to count the number of matches

  BlastMatch() {reset();};

  BlastMatch(const BlastMatch& align)
    {                   
      query_num=align.query_num;          
      subject_num=align.subject_num;      
      query_strand=align.query_strand;
      subject_strand=align.subject_strand;
      score=align.score;
      e_value=align.e_value;
      identity=align.identity;
      length=align.length;
      id=align.id;
      query_strt=align.query_strt;
      subject_strt=align.subject_strt;
      query_end=align.query_end;
      subject_end=align.subject_end;
    };
  
  
  virtual ~BlastMatch(){};
  
  virtual void* clone() const
    { 
      return (void*)new BlastMatch(*this);
    };
  
  unsigned long getQuery_num()              const   {return query_num;};
  unsigned long getSubject_num()            const   {return subject_num;};
  char getQuery_strand()           const   {return query_strand;};
  char getSubject_strand()         const   {return subject_strand;};
  unsigned long getScore()         const   {return score;};
  double getE_value()              const   {return e_value;};
  double getIdentity()             const   {return identity;};
  unsigned getLength()             const   {return length;};
  unsigned long getId()            const   {return id;};
  unsigned long getQuery_strt()    const   {return query_strt;};
  unsigned long getSubject_strt()  const   {return subject_strt;};
  unsigned long getQuery_end()     const   {return query_end;};
  unsigned long getSubject_end()   const   {return subject_end;};
  
  
  void setQuery_strand(char qd)           {query_strand=qd;};
  void setSubject_strand(char sd)         {subject_strand=sd;};
  void setQuery_num(unsigned long qn)              {query_num=qn;};
  void setSubject_num(unsigned long sn)            {subject_num=sn;};
  void setScore(unsigned long s)          {score=s;};
  void setE_value(double pv)              {e_value=pv;};
  void setIdentity(double i)              {identity=i;};
  void setLength(unsigned l)              {length=l;};
  void setId(unsigned long i)             {id=i;};
  void setQuery_strt(unsigned long  qs)   {query_strt=qs;};
  void setSubject_strt(unsigned long ss)  {subject_strt=ss;};
  void setQuery_end(unsigned long qe)     {query_end=qe;};
  void setSubject_end(unsigned long se)   {subject_end=se;};
  
  
  void reset(void);

  void readncbiblastfield(std::istream&);
  void readwublastfield(std::istream&);

  // Fonction to read binary file list of blast matches
  void readlst(std::istream&);
  void readoldlst(std::istream&);

  // Fonction to write binary file list of blast matches
  void writelst(std::ostream&);
  void view(void) const;
  // overload operateur >>
  friend std::istream& operator>>(std::istream&, BlastMatch&);
  // overload operateur >>
  friend std::istringstream& operator>>(std::istringstream&, BlastMatch&);
  // overload operateur <<
  friend std::ostream& operator<<(std::ostream&, BlastMatch&);
 
}; 

#endif









