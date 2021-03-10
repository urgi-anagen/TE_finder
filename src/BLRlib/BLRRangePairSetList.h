/*
 * BLRRangePairSetList.h
 *
 *  Created on: 14 oct. 2015
 *      Author: hquesnev
 */

#ifndef BLRRANGEPAIRSETLIST_H_
#define BLRRANGEPAIRSETLIST_H_

#include "RangePairSet.h"

class  BLRRangePairSetList : public std::list<RangePairSet>
{

public:
    BLRRangePairSetList(const std::list<RangePairSet>& rps):std::list<RangePairSet>(rps){};
    void writePath(std::ostream& out);
    void writePathAttr(std::ostream& out);
    void writeGFF3(std::ostream& out);

};


#endif /* BLRRANGEPAIRSETLIST_H_ */
