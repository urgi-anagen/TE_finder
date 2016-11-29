/*
 * BLRRangePairSetListWriter.cpp
 *
 *  Created on: 14 oct. 2015
 *      Author: hquesnev
 */
#include <stdlib.h>
#include <fstream>
#include "SDGError.h"

#include "BLRRangePairSetListWriter.h"

BLRRangePairSetListWriter::BLRRangePairSetListWriter()
{
}

BLRRangePairSetListWriter::~BLRRangePairSetListWriter()
{
}


//---------------------------------------------------------------------------------------
void BLRRangePairSetListWriter::writePath(const BLRRangePairSetList& rpsList,
		const SDGBioSeqDB& query_db, const SDGBioSeqDB& subject_db,
		std::ostream& out)
{
	unsigned path_id=0;
	for(BLRRangePairSetList::const_iterator iter_list
			=rpsList.begin();iter_list!=rpsList.end();
			iter_list++)
		{
		  unsigned id = ++path_id;
		  iter_list->write(out,id,query_db[iter_list->getNumQuery()-1].getDE(),subject_db[iter_list->getNumSubject()-1].getDE());
		}
}
//---------------------------------------------------------------------------------------
void BLRRangePairSetListWriter::writeGFF3(const BLRRangePairSetList& rpsList,
		const SDGBioSeqDB& query_db, const SDGBioSeqDB& subject_db,
		std::ostream& out, const std::string& source)
{
	unsigned path_id=0;
	for(BLRRangePairSetList::const_iterator iter_list
			=rpsList.begin();iter_list!=rpsList.end();
			iter_list++)
		{
		  unsigned id = ++path_id;
		  iter_list->writeGFF3(out,id,query_db[iter_list->getNumQuery()-1].getDE(),subject_db[iter_list->getNumSubject()-1].getDE(),source);
		}
}

