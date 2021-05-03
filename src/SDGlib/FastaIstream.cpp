#include <FastaIstream.h>
#include <BioSeq.h>
#include <string>
#include <sstream>

FastaIstream& FastaIstream::operator>>(BioSeq &p)
{
    std::string titre;
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
                        titre=line.substr(1);
                        titre=titre.erase(titre.find_last_not_of("\t\n\v\f\r ") + 1);
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

  p=BioSeq(buffstr.str(),titre);
  return *this;
}

