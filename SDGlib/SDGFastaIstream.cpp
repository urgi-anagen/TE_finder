#include <SDGFastaIstream.h>
#include <SDGString.h>
#include <sstream>


SDGFastaIstream& SDGFastaIstream::operator>>(SDGBioSeq &p)
{
    SDGString titre;
    std::string line;
    std::ostringstream buffstr;
    bool found=false;
    if (this->is_open())
      {
        while ( *this )
        {
            if(peek()=='>')
            {
                if(!found)
                {
                    if(std::getline(*this,line))
                    {
                        titre=SDGString(line).substr(1);
                        titre=titre.trimR();
                    }
                    found=true;
                }
                else
                    break;
            }
            else
            {
                if(std::getline(*this,line))
                    buffstr<<line;
            }          
          }         
        }

  p=newSDGMemBioSeq(buffstr.str());
  p.setDE(titre);

  return *this;
}

SDGFastaIstream SDGFastaCin(0);
