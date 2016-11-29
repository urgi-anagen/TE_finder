
#include <SDGBioSeqDB.h>
#include <SDGFastaBioSeq.h>
#include <fstream>


//-------------------------------------------------------
void SDGBioSeqDB::load( SDGString fichier, int verbose )
{
  clear();
  filename=fichier;
  if(verbose>0)
	  std::cout<<"Loading file '"<<fichier<<"'... "<<std::flush;
  std::ifstream file(fichier);
  unsigned count=0;
  while (file)
    {
      if(file.peek()=='>')
	{
	  push_back(newSDGFastaBioSeq(fichier,file.tellg()));
	  name2pos[back().getDE()]=count++;
	}
      file.get();
    }
  if(verbose>0)
	  std::cout<<"File was loaded."<<std::endl;

  file.close();
}
