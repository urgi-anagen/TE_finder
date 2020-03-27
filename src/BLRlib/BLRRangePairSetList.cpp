/*
 * BLRRangePairSetList.cpp
 *
 *  Created on: 14 oct. 2015
 *      Author: hquesnev
 */
#include <stdlib.h>
#include <fstream>
#include "SDGError.h"

#include "BLRRangePairSetList.h"



//---------------------------------------------------------------------------------------
void BLRRangePairSetList::writePath(std::ostream& out)
{
	unsigned path_id=0;
	for(BLRRangePairSetList::const_iterator iter_list
			=begin();iter_list!=end();
			iter_list++)
		{
		  unsigned id = ++path_id;
		  iter_list->write(out,id,iter_list->getRangeQ().getNameSeq(),iter_list->getRangeS().getNameSeq());
		}
}
//---------------------------------------------------------------------------------------
void BLRRangePairSetList::writePathAttr(std::ostream& out)
{
    unsigned path_id=0;
    for(BLRRangePairSetList::const_iterator iter_list
            =begin();iter_list!=end();
        iter_list++)
    {
        unsigned id = ++path_id;
        iter_list->writeRpsAttr(out,id,iter_list->getRangeQ().getNameSeq(),iter_list->getRangeS().getNameSeq());
    }
}
//---------------------------------------------------------------------------------------
void BLRRangePairSetList::writeGFF3(std::ostream& out, const std::string& source)
{
	unsigned path_id=0;
	for(BLRRangePairSetList::const_iterator iter_list
			=begin();iter_list!=end();
			iter_list++)
		{
		  unsigned id = ++path_id;
		  iter_list->writeGFF3(out,id,iter_list->getRangeQ().getNameSeq(),iter_list->getRangeS().getNameSeq(),source);
		}
}

