/*
 * Test_FragJoin.h
 *
 *  Created on: 11 oct. 2021
 *      Author: hquesnev
 */

#ifndef TEST_FRAGJOIN_H_
#define TEST_FRAGJOIN_H_
#include <cppunit/extensions/HelperMacros.h>

#include "FragJoin.h"
#include <list>
#include <vector>

class Test_FragJoin : public CppUnit::TestFixture {

	CPPUNIT_TEST_SUITE(Test_FragJoin);

	CPPUNIT_TEST( test_align_all );

	CPPUNIT_TEST_SUITE_END();

public:

    Test_FragJoin(void) {}

	void setUp()
	{
	}
	void tearDown()
	{
	}

protected:
    void test_align_all(void);

};


#endif /* TEST_FRAGJOIN_H_ */
