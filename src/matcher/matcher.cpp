/***
 * 
 * matcher.cpp
 *
 ***/

#include <iostream>

#include "BLRMatchMap.h"
#include "SDGString.h"
#include "Reference.h"

int main(int argc, char* argv[])
{
   try{
     std::cout<<"Beginning Matcher (version "<<VERSION<<")"<<std::endl<<std::flush;

     BLRMatcherThreadsParameter para;
     para.parseOptArg(argc,argv);
     if(para.getVerbose()>0)
       para.view(std::cout);

     BLRMatchMap match_map(para);

     if(para.getVerbose()>0)
       std::cout<<"Load the matches..."<<std::endl<<std::flush;
     match_map.load( para.getVerbose()-1 );
     if(para.getVerbose()>0)
       std::cout<<"Matches were loaded."<<std::endl;

     if(para.getCleanBefore())
       {
	 if(para.getVerbose()>0)
	   std::cout<<"Clean the conflicts..."<<std::endl<<std::flush;
	 match_map.clean_conflicts();
	 if(para.getVerbose()>0)
	   std::cout<<"Conflicts were cleaned."<<std::endl;
       }

     if(para.getJoin_frag())
       {
	 if(para.getVerbose()>0)
	   std::cout<<"Connect the fragments..."<<std::endl<<std::flush;
	 match_map.mapPath(true,para.getCleanBefore(),para.getCleanAfter(), para.getMerge(), para.getVerbose()-1);
	 if(para.getVerbose()>0)
	   std::cout<<"Fragments were connected."<<std::endl;
     
       }
     else
       match_map.mapPath(false);

     std::ostringstream filename;
     if(para.getCleanBefore() || para.getCleanAfter() )
       filename<<para.getPrefixFileName()<<".clean_match";
     else
       filename<<para.getPrefixFileName()<<".match";

     if(para.getVerbose()>0)
       std::cout<<"Write the results..."<<std::endl<<std::flush;

     RangeMap matchmap=match_map.writeMap(filename.str()+".map",para.getVerbose()-1);
     
     std::list<RangePairSet> rps_list = match_map.getRpsListFromMapPath();
     SDGString parafile=filename.str()+".param";
     para.write(filename.str()+".param");
     match_map.writePath(filename.str()+".path", rps_list, para.getVerbose()-1);
     match_map.writeBED(filename.str()+".bed", rps_list, para.getVerbose()-1);
     match_map.writeGFF3(filename.str()+".gff3", rps_list, para.getVerbose()-1);
     
     if(para.getBank()!="<not set!>" && para.getQuery()!="<not set>")
       match_map.writeMatch(filename.str()+".tab",para.getVerbose()-1);
     if(para.getQuery()!="<not set>")
       match_map.writeSeq(matchmap,filename.str()+".fa",para.getVerbose()-1);
     if(para.getVerbose()>0)
       std::cout<<"Results were written."<<std::endl<<std::flush;

     std::cout<<"End Matcher (version "<<VERSION<<")"<<std::endl;
     exit(EXIT_SUCCESS);
   }// end try

  catch(SDGException e)
    {
      std::cout<<std::endl<<e.message;
    }
  catch(...)
    {
      std::cout<<std::endl<<"unknown exception catch";
      exit(EXIT_FAILURE);
    }
}
