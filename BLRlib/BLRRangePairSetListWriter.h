/*
 * BLRRangePairSetListWriter.h
 *
 *  Created on: 14 oct. 2015
 *      Author: hquesnev
 */

#ifndef BLRRANGEPAIRSETLISTWRITER_H_
#define BLRRANGEPAIRSETLISTWRITER_H_

#include "SDGBioSeqDB.h"
#include "BLRRangePairSetList.h"

class BLRRangePairSetListWriter {
public:
	BLRRangePairSetListWriter();
	virtual ~BLRRangePairSetListWriter();

	void writePath(const BLRRangePairSetList& rpsList,
			const SDGBioSeqDB& query_db, const SDGBioSeqDB& subject_db,
			std::ostream& out);

	void writeGFF3(const BLRRangePairSetList& rpsList,
			const SDGBioSeqDB& query_db, const SDGBioSeqDB& subject_db,
			std::ostream& out, const std::string& source);
};

#endif /* BLRRANGEPAIRSETLISTWRITER_H_ */
