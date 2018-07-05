/*
 * Test_SDGlib.h
 *
 *  Created on: 5 juil 2018
 *      Author: hquesnev
 */

#ifndef TEST_SDGLIB_H_
#define TEST_SDGLIB_H_
#include <cppunit/extensions/HelperMacros.h>

#include <list>
#include <vector>

class Test_SDGlib : public CppUnit::TestFixture {

	CPPUNIT_TEST_SUITE(Test_SDGlib);
	
	CPPUNIT_TEST( test_SDGMemBioSeq );
	CPPUNIT_TEST( test_SDGFastaBioSeq );
	CPPUNIT_TEST( test_SDGFastaIstream );

	CPPUNIT_TEST_SUITE_END();

public:

	Test_SDGlib(void) {}

	void setUp()
	{}
	void tearDown()
	{}

protected:
	void test_SDGMemBioSeq(void);
	void test_SDGFastaBioSeq( void );
	void test_SDGFastaIstream( void );

};


#endif /* TEST_SDGLIB_H_ */
