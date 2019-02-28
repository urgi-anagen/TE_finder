#include "RangeMap.h"

#include <SDGError.h>
#include <SDGString.h>


int main(int argc, char* argv[])
{
  try{

    RangeMap range_container;

    if(argc!=4)
      {
	std::cerr<<"usage:"<<argv[0]<<" <map file> <fasta db> <size>"<<std::endl;
	exit(1);
      }
    unsigned size=atoi(argv[3]);
    std::cout<<"load map ..."<<std::flush;
    range_container.load(argv[1]);
    std::cout<<"done!"<<std::endl<<std::flush;
    range_container.view();
    unsigned countRange=range_container.getCountRange();
    std::cout<<"\t"<<countRange<<" found"<<std::endl;

    SDGString dbname(argv[1]);
    dbname=dbname.afterlast("/")+".flank5.fa";
    SDGBioSeqDB db(argv[2]);
    range_container.writeFlank5Seq(dbname,db,size);
  }
  catch(SDGException e)
    {
      std::cerr<<e.message<<std::endl;
    }
  return 1;
};



