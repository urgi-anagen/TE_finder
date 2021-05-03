/*
 * Test_BLRGroup.h
 *
 *  Created on: 30 déc. 2009
 *      Author: hquesnev
 */

#ifndef TEST_BLRGROUP_H_
#define TEST_BLRGROUP_H_
#include <cppunit/extensions/HelperMacros.h>
#include "Range.h"
#include "RangeAlign.h"
#include "RangeAlignSet.h"
#include "../BLRGroupThreads.h"
#include <list>
#include <vector>

class Test_BLRGroup : public CppUnit::TestFixture {

//	BLRGroup blrgroup;
//	BLRGrouperParameter para;
//	BLRMatchMap match_map;

	CPPUNIT_TEST_SUITE(Test_BLRGroup);
	
	CPPUNIT_TEST( test_eraseMember );
//	CPPUNIT_TEST( test_mergeWithExtGroup );

	CPPUNIT_TEST_SUITE_END();

public:

	Test_BLRGroup(void) {};

	void setUp()
	{
	}
	void tearDown()
	{
	}

protected:

	void test_eraseMember( void );
//	void test_mergeWithExtGroup(void);


};


#endif /* TEST_BLRGROUP_H_ */
