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

	BLRGroup* blrgroup_ptr;
	BLRGrouperParameter* para_ptr;
	BLRMatchMap* match_map_ptr;

	CPPUNIT_TEST_SUITE(Test_BLRGroup);
	
	CPPUNIT_TEST( test_eraseMember );
	CPPUNIT_TEST( test_mergeWithExtGroup );

	CPPUNIT_TEST_SUITE_END();

public:

	Test_BLRGroup(void) : blrgroup_ptr(0), para_ptr(0), match_map_ptr(0){}

	void setUp()
	{
	    para_ptr=new BLRGrouperParameter();
	    match_map_ptr=new BLRMatchMap(para_ptr);
	    blrgroup_ptr=new BLRGroup(para_ptr,match_map_ptr,match_map_ptr->getNbQseq());
	}
	void tearDown()
	{
		delete blrgroup_ptr;
		delete match_map_ptr;
		delete para_ptr;
	}

protected:

	void test_eraseMember( void );
	void test_mergeWithExtGroup(void);

};


#endif /* TEST_BLRGROUP_H_ */
