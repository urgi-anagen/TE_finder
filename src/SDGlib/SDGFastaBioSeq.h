/*
 *
 * Declaration d'une classe permettant de stocker une sequence biologique 
 *
 *
 */


#ifndef SDGFASTABIOSEQ_H
#define SDGFASTABIOSEQ_H

#include <ios>
#include <iostream>
#include <SDGString.h>
#include <SDGMemBioSeq.h>

#include <SDGReference.h>

class __SDGFastaBioSeq :  public __SDGBioSeq
{ 
  private:

  SDGString litSeq(unsigned long from,unsigned long lg) const;

  SDGString fichier;
  char forme;
  SDGString buffseq;
  unsigned buffseq_size;
  unsigned long buffseq_pos;
  std::streampos offset_sequence;
  unsigned npl;             // Nucleotide par ligne
  unsigned long longueur;
  __SDGBioSeq::type_molecule molecule;

  SDGString    banque;
  SDGString 	access;
  SDGString	identificateur;
  SDGString 	definition;

  public:
  
  __SDGFastaBioSeq(SDGString fich, unsigned long decalage=0, 
		 __SDGBioSeq::type_molecule mol=__SDGBioSeq::NDF,
		   unsigned buffsize=1000);

  __SDGFastaBioSeq(const __SDGFastaBioSeq &seq);

  virtual ~__SDGFastaBioSeq(void){};

  long getOffset() const { return offset_sequence; };

  virtual SDGString name() const { return "SDGFastaBioSeq"; };



    //Retourne un objet __SDGString contenant la sequence sous forme ascii
  virtual SDGString toString(unsigned long d=0, unsigned long lg=0) const;

  virtual std::string::const_iterator begin(void) const
    { throw SDGException(NULL,"SDGFastaBioSeq: begin() not implemented");};
  virtual std::string::const_iterator end(void) const
    { throw SDGException(NULL,"SDGFastaBioSeq: end() not implemented");};
  virtual std::string::const_reverse_iterator rbegin(void) const
    { throw SDGException(NULL,"SDGFastaBioSeq: rbegin() not implemented");};
  virtual std::string::const_reverse_iterator rend(void) const
    { throw SDGException(NULL,"SDGFastaBioSeq: rend() not implemented");};


  virtual SDGBioSeq subseq(unsigned long d, unsigned long lg=0)  const;
  
  virtual SDGBioSeq complement()  const ;
  virtual SDGBioSeq reverse()  const ;
  
  virtual SDGBioSeq concat(SDGBioSeq seq) const;
  
  virtual SDGBioSeq concat(SDGString seq) const;
  
  virtual SDGBioSeq concat(char base) const;

  virtual char charAt(unsigned long indice) const;

  void write_idx(std::ostream&);
  void read_idx(std::istream&);

  virtual unsigned long length() const ;

  bool isLinear() const;
  void setLinear(bool);

  virtual void *clone() const ;

	// retourne des informations sur la sequence
  virtual SDGString getBQ(void)  const ;
  virtual SDGString getAC(void)   const ;
  virtual SDGString getID(void)   const ;
  virtual SDGString getDE(void)   const ;
  virtual type_molecule     getMO(void)  const  ;

	// Permet de modifier les champs privees definissant la sequence   
  virtual void setBQ( const SDGString new_BQ);
  virtual void setAC( const SDGString new_AC);
  virtual void setID( const SDGString new_ID);
  virtual void setDE( const SDGString new_DE);

};


inline bool  __SDGFastaBioSeq::isLinear() const
{
  return forme;
};
inline void  __SDGFastaBioSeq::setLinear(bool x)
{
  forme = x;
};

inline SDGBioSeq newSDGFastaBioSeq(SDGString fich, long decalage=0, 
			    __SDGBioSeq::type_molecule mol=__SDGBioSeq::NDF)
{
  return SDGBioSeq(__SDGFastaBioSeq(fich,decalage,mol));
}

      
#endif










