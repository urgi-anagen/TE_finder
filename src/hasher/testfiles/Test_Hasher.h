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

class Test_Hasher : public CppUnit::TestFixture {

	CPPUNIT_TEST_SUITE(Test_Hasher);

    CPPUNIT_TEST( test_runAsScript );

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

    void test_runAsScript( void );

};


#endif /* TEST_DUSTER_H_ */
