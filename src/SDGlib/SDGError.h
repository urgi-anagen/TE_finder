#ifndef SDGERROR_H
#define SDGERROR_H

#include <exception>
#include <iostream>
#include <string>
#include <sstream>

void SDGError(int code);

class SDGException 
{
  public:

  const void *object;
  std::string message;
  int  error;

  SDGException(const void *o=NULL, std::string m="", int e=0) :
    object(o),
    message(m),
    error(e)
    {};
};

class Unsigned_Out_Of_Range : public std::exception{

private:
    std::string message;
    unsigned value;
    std::string what_string;

public:
    Unsigned_Out_Of_Range(std::string const& m="", unsigned v=0 ): message(m),value(v)
    {
        std::ostringstream sentence;
        sentence << message << value << " out of range !";
        what_string.assign(sentence.str());
    }

    virtual const char* what() const throw()
    {
        return what_string.c_str();
    }

};

class Long_Out_Of_Range : public std::exception{

private:
    std::string message;
    long value;
    std::string what_string;

public:
    Long_Out_Of_Range(std::string const& m="", long v=0 ): message(m),value(v)
    {
        std::ostringstream sentence;
        sentence << message << value << " out of range !";
        what_string.assign(sentence.str());
    }

    virtual const char* what() const throw()
    {
        return what_string.c_str();
    }

};
#endif /* SDGERROR_H */
