#ifndef TEST_ORIENTER_H
#define TEST_ORIENTER_H

#include <stdlib.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "SDGBioSeqDB.h"

class Test_orienter: public CPPUNIT_NS::TestFixture
{

	CPPUNIT_TEST_SUITE(Test_orienter);
        CPPUNIT_TEST(alignInputSequencesNotOnComplementWithoutScoreStorageAndSearchOnly);
        //CPPUNIT_TEST(hashAlign_search);
  	CPPUNIT_TEST_SUITE_END();

  	public:
		/*
     		void setUp();
		void tearDown();
		*/
	protected:
		void alignInputSequencesNotOnComplementWithoutScoreStorageAndSearchOnly(void);
		//void hashAlign_search(void);

	private:
		SDGBioSeqDB initBioSeqDB(void);
		std::vector<std::vector<int> > initVectorOfMatches(int lendb);
};
#endif
