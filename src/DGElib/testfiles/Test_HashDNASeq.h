/*
 * Test_Duster.h
 *
 *  Created on: 9 nov. 2015
 *      Author: hquesnev
 */

#ifndef TEST_HASHDNASEQ_H_
#define TEST_HASHDNASEQ_H_
#include <cppunit/extensions/HelperMacros.h>

#include "HashDNASeq.h"
#include <list>
#include <vector>

class Test_HashDNASeq : public CppUnit::TestFixture {

	CPPUNIT_TEST_SUITE(Test_HashDNASeq);
	
	CPPUNIT_TEST( test_hashSeqCount );
	CPPUNIT_TEST( test_hashSeqCountwHole );
	CPPUNIT_TEST( test_reverse_hash );
	CPPUNIT_TEST( test_reverse_hashwHole );
	CPPUNIT_TEST( test_diagSearch );

	CPPUNIT_TEST_SUITE_END();

public:

    Test_HashDNASeq(void) {}

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
};


#endif /* TEST_DUSTER_H_ */
