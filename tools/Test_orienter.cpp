#include "SDGString.h"
#include "SDGBioSeqDB.h"
#include "Test_orienter.h"
#include "HashAlignClone.h"

CPPUNIT_TEST_SUITE_REGISTRATION( Test_orienter );

SDGBioSeqDB Test_orienter::initBioSeqDB(void){
	SDGString filename = "./data/seqCluster2.fa";
	SDGBioSeqDB db( filename );
	return db;
}

void Test_orienter::alignInputSequencesNotOnComplementWithoutScoreStorageAndSearchOnly(void)
{

		HashAlignClone al;
		SDGBioSeqDB db = initBioSeqDB();
		bool verbose = false;
		for( SDGBioSeqDB::iterator i=db.begin(); i!=db.end();i++ )
		{
			for( SDGBioSeqDB::iterator j=i+1; j!=db.end();j++)
			{
  				bool b1 = i->getDE() == "BlastclustCluster2Mb22_chunk1617 (dbseq-nr 1) [174895,179677]";
  				bool b2 = j->getDE() == "BlastclustCluster2Mb35_chunk2038 (dbseq-nr 2) [84935,89713]";
				if (b1 && b2)
					verbose = true;

				al.setSeq(*i,*j);
				al.search();
			  	if (verbose)
					std::cout<<"E N D  F O R 1"<<std::endl<<std::flush;
 		
			}
			if (verbose)
				std::cout<<"E N D  F O R 2"<<std::endl<<std::flush;

		}
		if (verbose)
			std::cout<<"O U T  F O R 2"<<std::endl<<std::flush;


}

