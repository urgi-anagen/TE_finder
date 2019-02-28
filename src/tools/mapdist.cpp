#include <fstream>
#include <string>
#include <map>
#include <list>

#include <SDGString.h>
#include <RangeMap.h>

int main(int argc, char* argv[])
{
  try{

    RangeMap rangeMap;

    if(argc!=2)
      {
	std::cerr<<"usage:"<<argv[0]<<" <map file>"<<std::endl;
	exit(1);
      }


    std::cout<<"loading sequences range "<<argv[1]<< " ... "<<std::endl;
    rangeMap.load(argv[1]);
    std::cout<<"sequences range loaded !!"<<std::endl;

    rangeMap.sort();
    rangeMap.view();
  }
  catch(SDGException e)
    {
      std::cerr<<e.message<<std::endl;
    }
  return 1;
};


