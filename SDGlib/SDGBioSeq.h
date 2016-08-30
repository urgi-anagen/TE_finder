/*
 *
 * Declaration d'une classe permettant de stocker une sequence biologique 
 *
 *
 */


#ifndef SDGBIOSEQ_H
#define SDGBIOSEQ_H

#include <SDGString.h>
#include <SDGReference.h>
#include <SDGSequence.h>

class SDGBioSeq;

class __SDGBioSeq : public __SDGSequence
{ 

  protected :

	  /**
	   * Alphabet decrivant les proteines.
	   */
    static const std::string AAalphabet;

	  /**
	   * Alphabet decrivant l'ADN.
	   */
    static const std::string DNAalphabet;

	  /**
	   * Alphabet decrivant les sequences d'ADN degeneree.
	   */
    static const std::string DNADGalphabet;

    void checkPos(unsigned long &p) const;

 public:
  /**
   * enumeration des differents type de sequence gerable par une SDGBioSeq.
   * NDF = Type indefini <BR>
   * NUC = Acide nucleique <BR>
   * PRT = Proteine <BR>
   * NDG = Acide nucleique degenere
   * @see getMO
   * @see SDGBioSeq::getMO
   */
    enum type_molecule {NDF=0,NUC,PRT,NDG}; 
    static type_molecule sequenceType(const char *seq);

    static char complement(char a);	// ok

	virtual ~__SDGBioSeq() {};

	//Retourne un objet __SDGString contenant la sequence sous forme ascii
	virtual SDGString toString(unsigned long d=0,unsigned long lg=0) const =0;
	virtual std::string::const_iterator begin(void) const =0;
	virtual std::string::const_iterator end(void) const =0;
	virtual std::string::const_reverse_iterator rbegin(void) const =0;
	virtual std::string::const_reverse_iterator rend(void) const =0;

	virtual SDGBioSeq subseq(unsigned long d,unsigned long lg=0) const =0;

	virtual SDGBioSeq complement() const =0;

	virtual SDGBioSeq concat(SDGBioSeq seq) const =0;
	virtual SDGBioSeq concat(SDGString seq) const =0;
	virtual SDGBioSeq concat(char base) const =0;

	virtual SDGString name() const =0;

	virtual unsigned alphabetSize () const 
	  {
	    switch (getMO())
	      {
	      case NUC : return 4;
	      case NDG : return 16;
	      case PRT : return 22;
	      default  : return 0;
	      }
	  };
	virtual unsigned symbolLength () const
	  {
	    return 1;
	  }; 

	virtual void *clone() const =0;

	virtual bool isLinear() const;
	virtual void setLinear(bool);

	

	// retourne des informations sur la sequence
	virtual SDGString getBQ(void) const =0;
        virtual SDGString getAC(void) const  =0;
	virtual SDGString getID(void) const  =0;
	virtual SDGString getDE(void) const  =0;
	virtual type_molecule     getMO(void) const  =0;

	// Permet de modifier les champs privees definissant la sequence   
	virtual void setBQ(const SDGString new_BQ)=0;
        virtual void setAC(const SDGString new_AC)=0;
	virtual void setID(const SDGString new_ID)=0;
	virtual void setDE(const SDGString new_DE)=0;

	// Permet d'acceder aux annotation d'une sequence


	virtual SDGBioSeq fivePrime(unsigned long pos, unsigned long lg) const;
	virtual SDGBioSeq threePrime(unsigned long pos, unsigned long lg) const;
	virtual unsigned long posInRef(unsigned long indice) const;
	virtual SDGBioSeq getBioseq() const;
	virtual unsigned long posInAbsoluteRef(unsigned long indice) const;
	virtual SDGBioSeq getAbsoluteBioseq() const;

	std::string getAlphabet() const
	  {
	    switch (getMO())
	      {
	      case NUC : return DNADGalphabet;
	      case NDG : return DNAalphabet;
	      case PRT : return AAalphabet;
	      default  : return "";
	      }	    
	  }
};


      
#endif






