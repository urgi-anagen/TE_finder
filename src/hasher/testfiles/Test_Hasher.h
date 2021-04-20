/*
 * Test_Duster.h
 *
 *  Created on: 9 nov. 2015
 *      Author: hquesnev
 */

#ifndef TEST_HASHER_H_
#define TEST_HASHER_H_
#include <cppunit/extensions/HelperMacros.h>

#include "Hasher.h"
#include <list>
#include <vector>
#include <iomanip>

class Test_Hasher : public CppUnit::TestFixture {

	CPPUNIT_TEST_SUITE(Test_Hasher);

	CPPUNIT_TEST( test_search );
	CPPUNIT_TEST( test_diagSearchDist );
	CPPUNIT_TEST( test_diagSearchScore );

//	CPPUNIT_TEST( test_fragMerge );
//    CPPUNIT_TEST( test_runAsScript );

	CPPUNIT_TEST_SUITE_END();

public:

    Test_Hasher(void) {}

	void setUp()
	{
	}
	void tearDown()
	{
	}

protected:
    void test_search( void );
    void test_diagSearchDist( void );
    void test_diagSearchScore( void );
//    void test_fragMerge( void );
    void test_runAsScript( void );

};


#endif /* TEST_DUSTER_H_ */
