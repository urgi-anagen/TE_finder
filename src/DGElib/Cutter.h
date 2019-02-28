/**
 * \file Cutter.h
 * \brief Header file for the class Cutter
 */

#ifndef CUTTER_H
#define CUTTER_H

#include <stdlib.h>
#include "SDGString.h"
#include "Reference.h"

class Cutter
{
 private:
  unsigned word_default;        //!< word length default of blast (11)
  SDGString extention_default;  //!< extension file default (_cut)
  unsigned long length;         //!< Maximum length of subsequence
  unsigned over;                //!< Length of the overlaps between subsequence
  unsigned word;               /*!< Minimum length threshold of subsequence
			and Maximum length threshold of unknown subsequence*/
  SDGString extention;               //!< Extension for cut-out file
  int verbose;

 public:
  Cutter(){reset();};

  void reset(void)
    {
      word_default=11;
      extention_default="_cut";
      length=50000;
      over=100;
      word=word_default;
      extention=extention_default;
      verbose=0;
    };

  //! overload operator ==
  friend bool operator==(const Cutter&,const Cutter&);
  //! overload operator !=
  friend bool operator!=(const Cutter&,const Cutter&);

  /*! Function to get length of cut-out \return length*/
  unsigned long getLength() const                    {return length;};

  /*! Function to get over of cut-out \return over*/
  unsigned getOver()  const                     {return over;};

  /*! Function to get word of cut-out \return word*/
  unsigned getWord() const                      {return word;};

  /*! Function to get word_default value \return word_default*/
  unsigned getWord_default() const              {return word_default;};

  /*! Function to get extension of cut-out \return extension*/
  SDGString getExtention() const           {return extention;};

  /*! Function to get extention_default value \return extention_default*/
  SDGString getExtention_default() const   {return extention_default;};

  //! Function to set length
  void setLength(unsigned long l)              {length=l;};

  //! Function to set over
  void setOver(unsigned o)                {over=o;};

  //! Function to set word
  void setWord(unsigned w)                {word=w;};

  //! Function to set extension
  void setExtention(SDGString e)     {extention=e;};

  bool check( SDGString bank_name, int verbose=0 );

  /*! Function to cut databank*/
  SDGString cutDB( SDGString, int verbose=0 );
};

#endif
