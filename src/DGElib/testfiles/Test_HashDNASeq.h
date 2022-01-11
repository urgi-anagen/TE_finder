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
#include "MinimizerFuncDNASeq.h"
#include <list>
#include <vector>

class Test_HashDNASeq : public CppUnit::TestFixture {

	CPPUNIT_TEST_SUITE(Test_HashDNASeq);
	
	CPPUNIT_TEST( test_hashSeqCount );
	CPPUNIT_TEST( test_hashSeqCountwHole );
	CPPUNIT_TEST( test_reverse_hash );
	CPPUNIT_TEST( test_hash );
	CPPUNIT_TEST( test_reverse_hashwHole );
	CPPUNIT_TEST( test_diagSearchDist );
    CPPUNIT_TEST( test_minimizer );

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
    void test_hash(void);
	void test_reverse_hashwHole(void);
	void test_diagSearchDist( void );
    void test_minimizer(void);
};


#endif /* TEST_DUSTER_H_ */
