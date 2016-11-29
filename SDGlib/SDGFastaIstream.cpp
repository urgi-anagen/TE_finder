#include <SDGFastaIstream.h>
#include <SDGString.h>
#include <sstream>


SDGFastaIstream& SDGFastaIstream::operator>>(SDGBioSeq &p)
{
  char buff[1024];
  do{
    getline(buff,1024);
  }while((*this) && *buff!='>');
  buff[1023]='\0';
  SDGString titre=SDGString(buff).substr(1);
  titre=titre.trimR();
  std::ostringstream buffstr;
  while ( (*this) && peek()!='>')
    { 
      getline(buff,1024);
      buffstr<<buff;
    }
  p=newSDGMemBioSeq(buffstr.str());
  p.setDE(titre);

  return *this;
}

SDGFastaIstream SDGFastaCin(0);




