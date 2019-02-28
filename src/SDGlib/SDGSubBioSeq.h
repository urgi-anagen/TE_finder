/*
 *
 * Declaration d'une classe permettant de stocker une sequence biologique 
 *
 *
 */


#ifndef SDGSUBBIOSEQ_H
#define SDGSUBBIOSEQ_H

#include <vector>

#include <SDGString.h>
#include <SDGMemBioSeq.h>
#include <SDGReference.h>

class __SDGSubBioSeq : public __SDGBioSeq 
{ 
  private:

  SDGBioSeq bioseq;

  SDGString     banque;
  SDGString 	access;
  SDGString	identificateur;
  SDGString 	definition;
  __SDGBioSeq::type_molecule molecule;

 
  std::vector<unsigned long> debut, longueur;
  

  public:

  __SDGSubBioSeq(SDGBioSeq &seq, unsigned long d, 
		  unsigned long lg=0);
  
  __SDGSubBioSeq(const __SDGBioSeq &seq, unsigned long d, 
		  unsigned long lg=0);
  
  __SDGSubBioSeq(const __SDGSubBioSeq &seq);

  virtual ~__SDGSubBioSeq(void){};

  virtual SDGString name() const { return "SDGSubBioSeq"; }; 
  
  virtual SDGString toString(unsigned long d=0, unsigned long lg=0) const;

  virtual std::string::const_iterator begin(void) const
    { throw SDGException(NULL,"SDGSubBioSeq: begin() not implemented");};
  virtual std::string::const_iterator end(void) const
    { throw SDGException(NULL,"SDGSubBioSeq: end() not implemented");};

  virtual std::string::const_reverse_iterator rbegin(void) const
    { throw SDGException(NULL,"SDGSubBioSeq: rbegin() not implemented");};
  virtual std::string::const_reverse_iterator rend(void) const
    { throw SDGException(NULL,"SDGSubBioSeq: rend() not implemented");};


  virtual SDGBioSeq subseq(unsigned long d, unsigned long lg=0) const;
  
  virtual SDGBioSeq complement() const;
  
  virtual SDGBioSeq concat(SDGBioSeq seq) const;
  
  virtual SDGBioSeq concat(SDGString seq) const;
  
  virtual SDGBioSeq concat(char base) const;

  virtual char charAt(unsigned long indice) const;

  virtual unsigned long posInRef(unsigned long indice) const;

  virtual SDGBioSeq getBioseq() const;

  virtual unsigned long posInAbsoluteRef(unsigned long indice) const;

  virtual SDGBioSeq getAbsoluteBioseq() const;

  virtual long unsigned length() const;

  virtual void *clone() const ;

	// retourne des informations sur la sequence

  virtual SDGString getBQ(void) const;
  virtual SDGString getAC(void) const;
  virtual SDGString getID(void) const;
  virtual SDGString getDE(void) const;
  virtual __SDGBioSeq::type_molecule     getMO(void) const;

	// Permet de modifier les champs privees definissant la sequence   
  virtual void setBQ(const SDGString new_BQ);
  virtual void setAC(const SDGString new_AC);
  virtual void setID(const SDGString new_ID);
  virtual void setDE(const SDGString new_DE);


};

inline SDGBioSeq newSDGSubBioSeq(SDGBioSeq &seq, unsigned long d, 
		unsigned long lg=0)
{
  return SDGBioSeq(__SDGSubBioSeq(seq,d,lg));
};

inline SDGBioSeq newSDGSubBioSeq(const __SDGBioSeq &seq, unsigned long d, 
		unsigned long lg=0)
{
  return SDGBioSeq(__SDGSubBioSeq(seq,d,lg));
};


#endif








