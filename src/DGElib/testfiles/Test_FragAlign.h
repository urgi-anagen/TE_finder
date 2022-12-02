/*
 * Test_FragAlign.h
 *
 *  Created on: 11 oct. 2021
 *      Author: hquesnev
 */

#ifndef TEST_FRAGALIGN_H_
#define TEST_FRAGALIGN_H_
#include <cppunit/extensions/HelperMacros.h>

#include "FragAlign.h"
#include <list>
#include <vector>

class Test_FragAlign : public CppUnit::TestFixture {

	CPPUNIT_TEST_SUITE(Test_FragAlign);
	
	CPPUNIT_TEST( test_join );

	CPPUNIT_TEST_SUITE_END();

public:

    Test_FragAlign(void) {}

	void setUp()
	{
	}
	void tearDown()
	{
	}

protected:
	void test_join(void);

};


#endif /* TEST_FRAGALIGN_H_ */
