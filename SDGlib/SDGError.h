#ifndef SDGERROR_H
#define SDGERROR_H


#include <iostream>
#include <string>

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

#endif /* SDGERROR_H */
