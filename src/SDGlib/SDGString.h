/**
 * \file SDGString.cpp
 * \brief Header file for the class BLRGroup
 */

#ifndef SDGSTRING_H
#define SDGSTRING_H

#include <string>
#include <iostream>
#include <sstream>
#include <cctype>
#include <algorithm>
#include <vector>
#include <regex.h>

#include <SDGError.h>


/**
 * \class SDGString
 * \brief Implement a string with its alphabet
 */
class SDGString: public std::string
{
  protected :

  std::string alphabet;
  char defaultChar;
  void checkAlphabet(void)
    {
      if ( alphabet.empty()) return;

      for(std::string::iterator i=begin(); i!=end();i++)
    	  if(alphabet.find(*i)==std::string::npos)
    	  {
    		  std::ostringstream ostr;
    		  ostr<<*i<<" forbiden character!!";
    		  throw
			  	  SDGString::AlphabetException(this,ostr.str());
    	  }
    };


  public:

  /*******************************************
   *
   *  Constructeurs
   *
   *******************************************/

/*   SDGString(void):  */
/*     std::string(),alphabet(""),defaultChar('?') */
/*       {}; */
  SDGString(void):defaultChar('?')
      {};
  SDGString(const std::string& s, const std::string& alpha="", char def='?' ):
    std::string(s),alphabet(alpha),defaultChar(def)
    {};

  SDGString(const char* ch, const std::string& alpha="", char def='?' ):
    std::string(ch),alphabet(alpha),defaultChar(def)
    {};

  SDGString(char* ch, const std::string& alpha="", char def='?' ):
    std::string(),alphabet(alpha),defaultChar(def)
    {
      std::ostringstream ostr;
      ostr<<ch;
      *this=ostr.str();
    };

  SDGString(const char c, const std::string& alpha="", char def='?' ):
    std::string(),alphabet(alpha),defaultChar(def)
    {
      std::ostringstream ostr;
      ostr<<c;
      *this=ostr.str();
    };

  SDGString( double number, const std::string& alpha="", char def='?' ):
    std::string(),alphabet(alpha),defaultChar(def)
    {
      std::ostringstream ostr;
      ostr<<number;
      *this=ostr.str();
    };

  SDGString( long long  number, const std::string& alpha="", char def='?' ):
    std::string(),alphabet(alpha),defaultChar(def)
    {
      std::ostringstream ostr;
      ostr<<number;
      *this=ostr.str();
    };

  SDGString(unsigned long number, const std::string& alpha="", char def='?'):
    std::string(),alphabet(alpha),defaultChar(def)
    {
      std::ostringstream ostr;
      ostr<<number;
      *this=ostr.str();
    };

  SDGString(long number, const std::string& alpha="", char def='?' ):
    std::string(),alphabet(alpha),defaultChar(def)
    {
      std::ostringstream ostr;
      ostr<<number;
      *this=ostr.str();
    };

  SDGString(unsigned number, const std::string& alpha="", char def='?' ):
    std::string(),alphabet(alpha),defaultChar(def)
    {
      std::ostringstream ostr;
      ostr<<number;
      *this=ostr.str();
    };

  SDGString(int number, const std::string& alpha="", char def='?' ):
    std::string(),alphabet(alpha),defaultChar(def)
    {
      std::ostringstream ostr;
      ostr<<number;
      *this=ostr.str();
    };

  SDGString(short number, const std::string& alpha="", char def='?' ):
    std::string(),alphabet(alpha),defaultChar(def)
    {
      std::ostringstream ostr;
      ostr<<number;
      *this=ostr.str();
    };

  SDGString(unsigned short number, const std::string& alpha="", char def='?' ):
    std::string(),alphabet(alpha),defaultChar(def)
    {
      std::ostringstream ostr;
      ostr<<number;
      *this=ostr.str();
    };

/*   void* clone() const */
/*       {  */
/* 	return (void*)new SDGString(*this); */
/*       }; */


    /****************************************
     *
     * Classes d'exceptions
     *
     ****************************************/
  class Exception : public SDGException
    {
      public :
	Exception(const void *o = NULL, std::string m="", int e=0) :
	SDGException(o,m,e) {};
    };

  class MemoryException : public Exception
    {
      public :
	MemoryException(const void *o = NULL, std::string m="", int e=0) :
	Exception(o,m,e) {};
    };

  class AlphabetException : public Exception
    {
      public :
	AlphabetException(const void *o = NULL, std::string m="", int e=0) :
	Exception(o,m,e) {};
    };

  class EmptyException : public Exception
    {
      public :
	EmptyException(const void *o = NULL, std::string m="", int e=0) :
	Exception(o,m,e) {};
    };

  class OutOfRangeException : public Exception
    {
      public :
	OutOfRangeException(const void *o = NULL, std::string m="", int e=0) :
	Exception(o,m,e) {};
    };

  const std::string& getAlphabet() const
    { return alphabet; };


  void setAlphabet(const std::string& alpha, char def='?')
    {
      alphabet=alpha;
      defaultChar = def;
    };

    /**
     * retourne le caractere par defaut
     * @return le caractere par defaut
     * @see setDefaultChar
     */

    char getDefaultChar() const
      {return defaultChar;};


    // retourne le caractere a la position
    // decrite par indice

    char charAt(unsigned long indice) const
      {
	return operator[](indice);
      };

    SDGString substr(unsigned long pos, unsigned long l=0) const
      {
	if(l)
	  return SDGString(std::string::substr(pos,l),alphabet,defaultChar);
	else
	  return SDGString(std::string::substr(pos),alphabet,defaultChar);
      };

    const char* start(void) const
      {
	return c_str();
      };

    SDGString chop(void) const
      {
	return substr(0,length()-1);
      };

    SDGString trimL(void) const
      {
	unsigned long i=0;
	while (isspace(charAt(i++)));
	return substr(i-1);
      };

    SDGString trimR(void) const
      {
	unsigned long i=length();
	while (isspace(charAt(--i)));
	return substr(0,i+1);
      };

    SDGString  operator+(const char* s) const
      {
	return std::operator+((std::string)*this,std::string(s));
      };

    SDGString  operator+(const SDGString s) const
      {
	return std::operator+((std::string)*this,std::string(s));
      };

    SDGString  operator+(const std::string& s) const
      {
	return std::operator+((std::string)*this,s);
      };

    SDGString  operator+(long long d) const
      {
	return std::operator+((std::string)*this,(std::string)SDGString((long long)d));
      };

    SDGString  operator+(unsigned long d) const
      {
	return std::operator+((std::string)*this,(std::string)SDGString((long long)d));
      };

    SDGString  operator+(double d) const
      {
	return std::operator+((std::string)*this,(std::string)SDGString(d));
      };

    void  operator+=(const char* s)
    {
      (*this)=operator+(s);
    };

//    void  operator+=(const char c)
//    {
//      (*this)=operator+(c);
//    };

    void  operator+=(const std::string& s)
    {
      (*this)=operator+(s);
    };

    void  operator+=(long long d)
    {
      (*this)=operator+(SDGString((long long)d));
    };

    void  operator+=(long d)
    {
      (*this)=operator+(SDGString((long long)d));
    };

    void  operator+=(unsigned long d)
    {
      (*this)=operator+(SDGString((long long)d));
    };

    void  operator+=(unsigned d)
    {
      (*this)=operator+(SDGString((long long)d));
    };

    void  operator+=(int d)
    {
      (*this)=operator+(SDGString((long long)d));
    };
    void  operator+=(short d)
    {
      (*this)=operator+(SDGString((long long)d));
    };
    void  operator+=(unsigned short d)
    {
      (*this)=operator+(SDGString((long long)d));
    };

    void operator+=(double d)
    {
      (*this)=operator+(SDGString(d));
    };

    operator const char*() const
      {
	return c_str();
      };


    //
	// Fonction de test de match avec une expression reguliere
	//

    long match(const SDGString& pattern, long& matchlen, unsigned offset=0) const
      {
	regex_t preg;
	size_t  nmatch=1;
	regmatch_t   pmatch[1];
	int eflags=0;

	regcomp(&preg, pattern, REG_EXTENDED | REG_NEWLINE);
	regexec(&preg,start()+offset,nmatch,pmatch,eflags);

	regfree(&preg);

	matchlen=pmatch[0].rm_eo-pmatch[0].rm_so;
	return pmatch[0].rm_so;
      };

    SDGString beforematch
      (const SDGString& pattern, unsigned  offset=0) const
      {
	long len=0;
	long l= match(pattern,len,offset);
	return substr(0,l+offset);
      };

    SDGString aftermatch
      (const SDGString& pattern, unsigned offset=0) const
      {
	long len=0;
	long l= match(pattern,len,offset);
	return substr(l+len+offset);
      };

    SDGString beforematch
      (const char* pattern, unsigned  offset=0) const
      {
	long len=0;
	long l= match(SDGString(pattern),len,offset);
	return substr(0,l+offset);
      };

    SDGString aftermatch
      (const char* pattern, unsigned offset=0) const
      {
	long len=0;
	long l= match(SDGString(pattern),len,offset);
	return substr(l+len+offset);
      };

    SDGString before(const SDGString& str) const
      {
	size_type pos=find(str);
	if(pos>size())
	  return *this;
	return substr(0,pos);
      }
    SDGString after(const SDGString& str) const
      {
	size_type pos=find(str);
	if(pos>size())
	  return *this;
	return substr(pos,size()-pos);
      }

    SDGString beforelast(const SDGString& str) const
      {
	size_type pos=rfind(str);
	if(pos>size())
	  return *this;
	return substr(0,pos);
      }
    SDGString afterlast(const SDGString& str) const
      {
	size_type pos=rfind(str);
	if(pos>size())
	  return *this;
	return substr(pos+1,size()-pos);
      }

//    SDGString toUpper(void) const
//      {
//	std::string  capital_s;
//	capital_s.resize(size());
//	std::transform(begin(), end(), capital_s.begin(), toupper);
//	return SDGString(capital_s,alphabet,defaultChar);
//      };
//
//    SDGString toLower(void) const
//    {
//    	std::string minus_s;
//    	minus_s.resize( size() );
//    	std::transform( begin(), end(), minus_s.begin(), tolower );
//    	return SDGString( minus_s, alphabet, defaultChar );
//    };

    static std::vector<SDGString> tokenize( const SDGString& str, const SDGString& delimiters="\t" );

};

#endif
