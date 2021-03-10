/*
 * Test_Duster.h
 *
 *  Created on: 9 nov. 2015
 *      Author: hquesnev
 */

#ifndef TEST_DUSTER_H_
#define TEST_DUSTER_H_
#include <cppunit/extensions/HelperMacros.h>

#include "Duster.h"
#include <list>
#include <vector>
#include <iomanip>

class Test_Duster : public CppUnit::TestFixture {

	CPPUNIT_TEST_SUITE(Test_Duster);
	CPPUNIT_TEST( test_search );
	CPPUNIT_TEST( test_fragMerge );
    CPPUNIT_TEST( test_runAsScript );

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
    void test_search( void );
	void test_fragMerge( void );
    void test_runAsScript( void );

};


#endif /* TEST_DUSTER_H_ */
