/**
 *
 * BLRMatchMapWriter.cpp
 *
 **/
#include <stdlib.h>
#include <fstream>
#include "BLRMatchMapWriter.h"
#include <SDGError.h>

//---------------------------------------------------------------------------------------
void BLRMatchMapWriter::writePath(const RpsList& rpsList,
		std::map<long,std::string> num2nameQ, std::map<long,std::string> num2nameS,std::ostream& out)
{
	unsigned path_id=0;
	for(RpsList::iterator iter_list
			=rpsList.begin();iter_list!=rpsList.end();
			iter_list++)
		{
		  std::string query_name=num2nameQ[iter_list->getNumQuery()];
		  std::string subject_name;

		  if(same_db)
			  subject_name=num2nameQ[iter_list->getNumSubject()];
		  else
			  subject_name=num2nameS[iter_list->getNumSubject()];
		  unsigned id = ++path_id;
		  iter_list->write(out,id,query_name,subject_name);
		}
}
//---------------------------------------------------------------------------------------
void BLRMatchMapWriter::writeGFF3(const RpsList& rpsList,
		std::map<long,std::string> num2nameQ, std::map<long,std::string> num2nameS,std::ostream& out, const std::string& source)
{
	unsigned path_id=0;
	for(RpsList::iterator iter_list
			=rpsList.begin();iter_list!=rpsList.end();
			iter_list++)
		{
		  std::string query_name=num2nameQ[iter_list->getNumQuery()];
		  std::string subject_name;

		  if(same_db)
			  subject_name=num2nameQ[iter_list->getNumSubject()];
		  else
			  subject_name=num2nameS[iter_list->getNumSubject()];
		  unsigned id = ++path_id;
		  iter_list->writeGFF3(out,id,query_name,subject_name);
		}
}
