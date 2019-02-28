/**
 * \file Test_BLRGroup.h
 * \brief Unitary tests for class BLRGroup
 */

#ifndef TEST_BLRGROUP_H
#define TEST_BLRGROUP_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

class Test_BLRGroup: public CPPUNIT_NS::TestFixture
{

	CPPUNIT_TEST_SUITE( Test_BLRGroup );
	//CPPUNIT_TEST( test_group );
	//CPPUNIT_TEST( test_group_load_path );
	//CPPUNIT_TEST( test_compare_join_for_debug );
	CPPUNIT_TEST( test_group_compare_intermediate_rpsList_between_2_loads );
 	CPPUNIT_TEST_SUITE_END();
	
	public:
		void setUp();
		void tearDown();
	
	protected:
		void test_group(void);
		void test_group_load_path(void);
		void test_compare_join_for_debug(void);
		void test_group_compare_intermediate_rpsList_between_2_loads(void);
};


#endif
