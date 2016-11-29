/*
 *
 * Declaration d'une classe permettant de stocker une sequence biologique 
 *
 *
 */


#ifndef SDGSEQUENCE_H
#define SDGSEQUENCE_H

#include <sstream>
#include <SDGString.h>
#include <SDGReference.h>

class SDGSequence;

class __SDGSequence
{ 

  public:

	virtual ~__SDGSequence() {};

	//Retourne un objet __SDGString contenant la sequence sous forme ascii
	virtual SDGString toString(unsigned long d=0,unsigned long lg=0) const =0;

	virtual SDGString name() const =0;

	virtual unsigned long length() const =0;

	virtual char charAt(unsigned long indice) const =0;	    

	virtual void *clone() const =0;

	virtual unsigned alphabetSize() const =0;

	virtual unsigned symbolLength() const =0;

};

class SDGSequence : public SDGReference<__SDGSequence>
{
 public:
	virtual ~SDGSequence() {};

	/**
	 * Retourne une representation de la sequence formatee 
         * sous la forme d'une SDGString.
	 *
	 * @param d long int autorisant la specification de la position
         *          a partir de laquel la convertion est realisee (valeur
         *          par defaut 0)
	 * @param lg long int specifiant sur quelle longueur a partir de
         *          d la convertion est realisee (valeur par defaut -1, soit
         *          la longueur permetant de convertir la totalite de la
         *          sequence
	 * @return  un objet de type SDGString representant la sequence
	 *          sous forme ASCII a partir de la position d sur une
	 *          longueur l
	 */
	virtual SDGString toString(unsigned long d=0,unsigned long lg=0) const
	  {
	    return getPointer()->toString(d,lg);
	  };

	/**
	 * Retourne le nom du type de l'objet.
	 * 
	 * @return un objet de type SDGString repr�sentant sous
         *         forme ascii le nom du type de l'objet
	 *         reference par la classe SDGSequence
	 * @see    SDGData:name()
	 */
	virtual SDGString name() const
	  {
	    return getPointer()->name();
	  };

	/**
	 * Retourne la longueur de la sequence
	 * @return  un long int representant le nombre de
	 *          lettre present dans la sequence.
	 */
	virtual unsigned long length() const
	  {
	    return getPointer()->length();
	  };

	/**
	 * Operateur d'acces a une lettre particuliaire de
	 * la sequence. La lettre est represente sous la forme
	 * d'un code numerique
	 * @param  indice  long int specifiant la position demandee.
	 * @return un long int representant la la lettre a la 
	 *         position indice sous la forme de son code 
	 *         numerique.
	 * @see    symbolAt
	 * @see    symbolTable
	 */
	char operator[](unsigned long indice) const
	  {
	    return getPointer()->charAt(indice);
	  };

	/**
	 * Retourne la taille de l'alphabet decrivant la sequence.
	 * @return un long int representant la taille de l'alphabet
	 *         decrivant la sequence
	 * @see    symbolTable
	 */
	unsigned long alphabetSize() const 
	  {
	    return getPointer()->alphabetSize();
	  };

	/**
	 * Retourne la longeur du symbole ASCII represantant
	 * une des lettres de l'alphabet utilise pour la sequence
	 * @return un long int representant la longueur du symbole
	 * @see    symbolTable
	 */
	unsigned long symbolLength() const 
	  {
	    return getPointer()->symbolLength();
	  };

	/**
	 * Retourne un tableau de SDGString correspondant a
	 * chacun des symboles decrivant les lettres de l'alphabet
	 * utilisee dans cette sequence. Il y a une relation direct
	 * entre la position dans cette table et le code numerique
	 * de la lettre. Si une lettre peut posseder deux orientations,
	 * celle-ci est codee par le signe de la lettre. Dans ce cas
	 * le symbole correspondant � la lettre est stocke dans la case
	 * correspondant a la valeur absolue du code.
	 * @return un objet du type SDGStringVector contenant l'ensemble
         *         des symboles decrivant l'alphabet
	 * @see    alphabetSize
	 * @see    symbolLength
	 * @see    letterAt
	 * @see    symbolAt
	 */

};
      
#endif





