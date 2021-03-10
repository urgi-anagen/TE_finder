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
		const std::map<long,std::string>& num2nameQ, const std::map<long,std::string>& num2nameS, std::ostream& out)
{
	unsigned path_id=0;
	for(RpsList::const_iterator iter_list
			=rpsList.begin();iter_list!=rpsList.end();
			iter_list++)
		{
		  std::string query_name=num2nameQ.at(iter_list->getNumQuery());
		  unsigned id = ++path_id;
		  iter_list->write(out,id,query_name,num2nameS);
		}
}
//---------------------------------------------------------------------------------------
void BLRMatchMapWriter::writeGFF3(const RpsList& rpsList,
		const std::map<long,std::string>& num2nameQ, const std::map<long,std::string>& num2nameS, std::ostream& out)
{
	unsigned path_id=0;
	for(RpsList::const_iterator iter_list
			=rpsList.begin();iter_list!=rpsList.end();
			iter_list++)
		{
		  std::string query_name=num2nameQ.at(iter_list->getNumQuery());
		  unsigned id = ++path_id;
		  iter_list->writeGFF3(out,id,query_name,num2nameS);
		}
}
