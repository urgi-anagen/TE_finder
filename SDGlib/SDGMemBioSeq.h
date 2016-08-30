/*
 * Class <code>__SDGMemBioSeq</code>
 *
 * La classe __SDGMemBioSeq derive de la class SDGData et __SDGBioSeq elle 
 * permet la manipulation des sequence biologique en memoire
 * dynamiquement.<BR>
 *
 * @author  Eric Coissac
 * @version 0.1, 05/05/99
 * @since   SDG1.0 *
 */


#ifndef SDGMEMBIOSEQ_H
#define SDGMEMBIOSEQ_H

#include <SDGReference.h>
#include <SDGString.h>
#include <SDGBioSeq.h>
#include <SDGFastaOstream.h>


#define PROTEINE PRT



class __SDGMemBioSeq :   public __SDGBioSeq
{ 

  friend class __SDGFastaBioSeq;

  private:

  type_molecule molecule;            // Type de la molecule (ADN,PRT...)
  bool  forme;                 // Vrai -> lineaire
                                     // faux -> circulaire
 
  SDGString     sequence;

  SDGString     banque;
  SDGString 	access;
  SDGString	    identificateur;
  SDGString 	definition;

  public:

	__SDGMemBioSeq(const __SDGMemBioSeq &ch);			     
	__SDGMemBioSeq(SDGString ch="", type_molecule mol=NDG);   

	virtual ~__SDGMemBioSeq(void){};

	virtual SDGString name() const { return "SDGMemBioSeq"; }; 

    virtual SDGString toString(unsigned long d=0, unsigned long lg=0) const;

	virtual std::string::const_iterator begin(void) const 
	  {return sequence.begin();};

	virtual std::string::const_iterator end(void) const 
	  {return sequence.end();};

	virtual std::string::const_reverse_iterator rbegin(void) const 
	  {return sequence.rbegin();};

	virtual std::string::const_reverse_iterator rend(void) const 
	  {return sequence.rend();};

	SDGBioSeq complement() const ;

	SDGBioSeq concat(SDGBioSeq seq) const;
	SDGBioSeq concat(SDGString seq) const;
	SDGBioSeq concat(char base) const;

	char charAt(unsigned long indice) const
	  { 
	    checkPos(indice);
	    return sequence.charAt(indice);
	  };	    

	unsigned long length() const;

	bool isLinear() const;
	void setLinear(bool);

	SDGBioSeq subseq(unsigned long d, unsigned long lg=0) const;

	void *clone() const;

	// retourne des informations sur la sequence
	SDGString         getBQ(void) const;
    SDGString         getAC(void) const;
	SDGString         getID(void) const;
	SDGString         getDE(void) const;
	type_molecule     getMO(void)  const;

	// Permet de modifier les champs privees definissant la sequence   
	void setBQ(const SDGString new_BQ);
    void setAC(const SDGString new_AC);
	void setID(const SDGString new_ID);
	void setDE(const SDGString new_DE);

};

//==========================================================================

class SDGBioSeq : public SDGReference<__SDGBioSeq>
{
 public:
 
  /*********************************
   *
   * Constructeurs
   *
   *********************************/

  SDGBioSeq(const __SDGBioSeq &obj=__SDGMemBioSeq());
  SDGBioSeq(const SDGBioSeq &orig);

  /****************************************
   *
   * Destructeur
   *
   ****************************************/


  /****************************************
   *
   * Classes d'exceptions
   *
   ****************************************/
    class Exception : public SDGException 
      { 
	public :
	Exception(const void *o = NULL,std::string m="", int e=0) :
	  SDGException(o,m,e) {};
      };

    /**
     *
     * La classe NotDefinedException est g�n�r�e pour les fonctions
     * virtuelle non redefinie par certaines structures de donn�es
     * h�ritant de la classe SDGBioSeq
     */

    class NotDefinedException : public Exception 
      { 
	public :
	NotDefinedException(const void *o = NULL,std::string m="", int e=0) :
	  Exception(o,m,e) {};
      };

    class OutOfRangeException : public Exception 
      { 
	public :
	OutOfRangeException(const void *o = NULL,std::string m="", int e=0) :
	  Exception(o,m,e) {};
      };

    class EmptyException : public Exception 
      { 
	public :
	EmptyException(const void *o = NULL,std::string m="", int e=0) :
	  Exception(o,m,e) {};
      };

    /**
     * name of the class type.
     * @return the name of the class type
     * @see    SDGData::name
     */

  SDGString name() const { return getPointer()->name(); }; 

  /**
   * Create from the SDGBioSeq object an SDGString representation of the 
   * sequence.
   * @param d start position from which the sequence is converted (default value is 0
   *           the begining of the sequence.
   * @param lg length of the converted sequence. the default value -1 indicate an convertion
   *            to the end of the sequence.
   * @return an SDGString repesentation of the sequence.
   * @see    SDGData::toString
   */


  SDGString toString(unsigned long d=0,unsigned long lg=0) const
    {
      return getPointer()->toString(d,lg);
    };

  std::string::const_iterator begin(void) const
    {
      return getPointer()->begin();
    };
  std::string::const_iterator end(void) const
    {
      return getPointer()->end();
    };
  std::string::const_reverse_iterator rbegin(void) const
    {
      return getPointer()->rbegin();
    };
  std::string::const_reverse_iterator rend(void) const
    {
      return getPointer()->rend();
    };

  /**
   * Create an SDGBioSeq object corresponding to a subsequence of the sequence.
   * @param d start position from which the sequence is converted (default value is 0
   *           the begining of the sequence.
   * @param lg length of the converted sequence. the default value -1 indicate an convertion
   *            to the end of the sequence.
   * @return an SDGBioSeq corresponding to the subsequence.
   * @see   __SDGMemBioSeq::subseq
   * @see   __SDGFastaBioSeq::subseq
   * @see   __SDGSubBioSeq::subseq
   * @see   isLinear
   * @see   setLinear
   */

  SDGBioSeq subseq(unsigned long d,unsigned long lg=0)   
    {
      return getPointer()->subseq(d,lg);
    };


  /**
   * return the length of the sequence.
   * @return a long integer corresponding to the length of the sequence.
   */

  unsigned long length() const
    {
      return getPointer()->length();
    };


  /**
   * return true if the object represente a linear sequence, and false in case of
   * a circular one.
   * @return a bool object <UL><LI><B>true</B> correspond to a linear sequence
   *                                 <LI><B>false</B> correspond to a circular sequence</UL>
   * @see setLinear
   */
  bool isLinear() const
    {
      return getPointer()->isLinear();
    };

  /**
   * change the linear statut of the sequence.
   * param x a bool object <UL><LI><B>true</B> correspond to a linear sequence
   *                                 <LI><B>false</B> correspond to a circular sequence</UL>
   * @see isLinear
   */
  void setLinear(bool x)
    {
      getMutablePointer()->setLinear(x);
    };

  /**
   * if the SDGBioSeq object correspond to a __SDGSubBioSeq fivePrime fonction return an 
   * SDGBioSeq object corresponding to the sequence in 5' of the subsequence
   * @param p start position of the returned sequence relative to the 5' position of
   *           the __SDGSubBioSeq
   * @param lg length of the returned sub sequence
   * @return an SDGBioSeq instance correponding to the subsequence describe by the two
   *         paramettre <B>p</B> and <B>lg</B>
   */

  SDGBioSeq fivePrime(unsigned long p, unsigned long lg) const
    {
      return getPointer()->fivePrime(p,lg);
    };

  /**
   * if the SDGBioSeq object correspond to a __SDGSubBioSeq threePrime fonction return an 
   * SDGBioSeq object corresponding to the sequence in 3' of the subsequence
   * @param p start position of the returned sequence relative to the 3' position of
   *           the __SDGSubBioSeq
   * @param lg length of the returned sub sequence
   * @return an SDGBioSeq instance correponding to the subsequence describe by the two
   *         paramettre <B>p</B> and <B>lg</B>
   */

  SDGBioSeq threePrime(unsigned long p, unsigned long lg) const
    {
      return getPointer()->threePrime(p,lg);
    };

  unsigned long posInRef(unsigned long indice) const
    {
      return getPointer()->posInRef(indice);
    };

  SDGBioSeq getBioseq() const
    {
      return getPointer()->getBioseq();
    };

  unsigned long posInAbsoluteRef(unsigned long indice) const
    {
      return getPointer()->posInAbsoluteRef(indice);
    };

  SDGBioSeq getAbsoluteBioseq() const
    {
      return getPointer()->getAbsoluteBioseq();    
    };

  SDGBioSeq complement() const
    {
      return  getPointer()->complement();
    };


  /**
   * return the char at the position indice of the sequence.
   * @param indice position of the char in the sequence
   * @return a char corresponding to the position indice in the sequence
   */
  char charAt(const unsigned long indice) const 
    {
      return getPointer()->charAt(indice);
    };	    

  // retourne des informations sur la sequence
  SDGString getBQ(void) const
    {
      return getPointer()->getBQ();
    };
      
  SDGString getAC(void)  const
    {
      return getPointer()->getAC();
    };

  SDGString getID(void)  const
    {
      return getPointer()->getID();
    };

  SDGString getDE(void) const
    {
      return getPointer()->getDE();
    };

  __SDGBioSeq::type_molecule getMO(void) const  
    {
      return getPointer()->getMO();
    };

  // Permet de modifier les champs privees definissant la sequence   

  void setBQ(SDGString new_BQ)
    {
      getMutablePointer()->setBQ(new_BQ);
    };

  void setAC(SDGString new_AC)
    {
      getMutablePointer()->setAC(new_AC);
    };

  void setID(SDGString new_ID)
    {
      getMutablePointer()->setID(new_ID);
    };

  void setDE(SDGString new_DE)
    {
      getMutablePointer()->setDE(new_DE);
    };


  SDGBioSeq operator + (SDGBioSeq seq) const
    {
      return getPointer()->concat(seq);
    };
  SDGBioSeq operator + (SDGString seq) const
    {
      return getPointer()->concat(seq);
    };
  SDGBioSeq operator + (char base) const
    {
      return getPointer()->concat(base);
    };

  SDGBioSeq &operator += (SDGBioSeq seq)
    {
      *this = getPointer()->concat(seq);
      return *this;
    };

  SDGBioSeq &operator += (SDGString seq)
    {
      *this = getPointer()->concat(seq);
      return *this;
    };

  SDGBioSeq &operator += (char base)
    {
      *this = getPointer()->concat(base);
      return *this;
    };

  std::string getAlphabet(void) const {return getPointer()->getAlphabet();};
};

//==========================================================================




inline SDGString __SDGMemBioSeq::getBQ()  const
{
  return banque;
}

inline SDGString __SDGMemBioSeq::getAC()  const
{
  return access;
}

inline SDGString __SDGMemBioSeq::getID()  const
{
  return identificateur;
}

inline SDGString __SDGMemBioSeq::getDE()  const
{
  return definition;
}

inline __SDGBioSeq::type_molecule __SDGMemBioSeq::getMO()  const
{
  return molecule;
}
  
inline void __SDGMemBioSeq::setBQ(const SDGString new_BQ)
{
  banque = new_BQ;
}

inline void __SDGMemBioSeq::setAC(const SDGString new_AC)
{
  access = new_AC;
}

inline void __SDGMemBioSeq::setID(const SDGString new_ID)
{
  identificateur = new_ID;
}
       
inline void __SDGMemBioSeq::setDE(const SDGString new_DE)
{
  definition = new_DE;
}
       	    
inline unsigned long __SDGMemBioSeq::length() const
{
  return sequence.length();
}

inline bool  __SDGMemBioSeq::isLinear() const
{
  return forme;
};
inline void  __SDGMemBioSeq::setLinear(bool x)
{
  forme = x;
};

inline SDGBioSeq __SDGMemBioSeq::subseq(unsigned  long d, 
				      unsigned long lg)  const
{
  __SDGMemBioSeq tmp(toString(d,lg),getMO());
  tmp.molecule = getMO();
  tmp.banque = getBQ();
  tmp.access = getAC();
  tmp.identificateur = getID();
  tmp.definition = getDE();

  return SDGBioSeq(tmp);
}


//==========================================================================

inline  SDGBioSeq::SDGBioSeq(const __SDGBioSeq &obj) : 
  SDGReference<__SDGBioSeq>(obj)  {};

inline  SDGBioSeq::SDGBioSeq(const SDGBioSeq &orig) : 
    SDGReference<__SDGBioSeq>(orig)  {};

inline SDGBioSeq newSDGMemBioSeq(SDGBioSeq ch)
{
  __SDGMemBioSeq seq(ch.toString(),ch.getMO());
  seq.setAC(ch.getAC());
  seq.setID(ch.getID());
  seq.setDE(ch.getDE());
  seq.setBQ(ch.getBQ());
		     
  return SDGBioSeq(seq);
}

inline SDGBioSeq newSDGMemBioSeq(SDGString ch, 
			      __SDGBioSeq::type_molecule mol=__SDGBioSeq::NDG)
{
  return SDGBioSeq(__SDGMemBioSeq(ch,mol));
}


#endif









