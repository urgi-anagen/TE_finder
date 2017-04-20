/*
 * Test_Duster.h
 *
 *  Created on: 9 nov. 2015
 *      Author: hquesnev
 */

#ifndef TEST_DUSTER_H_
#define TEST_DUSTER_H_
#include <cppunit/extensions/HelperMacros.h>

#include "../Duster.h"
#include <list>
#include <vector>

class Test_Duster : public CppUnit::TestFixture {

	CPPUNIT_TEST_SUITE(Test_Duster);
	
	CPPUNIT_TEST( test_hashSeqCount );
	CPPUNIT_TEST( test_hashSeqCountwHole );
	CPPUNIT_TEST( test_reverse_hash );
	CPPUNIT_TEST( test_reverse_hashwHole );
	CPPUNIT_TEST( test_diagSearch );
	CPPUNIT_TEST( test_fragMerge );

	CPPUNIT_TEST_SUITE_END();

public:

	Test_Duster(void) {}

	void setUp()
	{
	}
	void tearDown()
	{
	}

protected:
	void test_hashSeqCount(void);
	void test_hashSeqCountwHole(void);
	void test_reverse_hash(void);
	void test_reverse_hashwHole(void);
	void test_diagSearch( void );
	void test_fragMerge( void );

};


#endif /* TEST_DUSTER_H_ */
