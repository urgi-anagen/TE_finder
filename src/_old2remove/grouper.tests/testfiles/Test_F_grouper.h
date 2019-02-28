/**
 * \file Test_F_grouper.h
 * \brief Functional tests for the GROUPER program
 */

#ifndef TEST_F_GROUPER_H
#define TEST_F_GROUPER_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "SDGString.h"

/**
 * \class Test_F_grouper
 * \brief Functional tests for the GROUPER program
 */
class Test_F_grouper: public CPPUNIT_NS::TestFixture
{
	CPPUNIT_TEST_SUITE( Test_F_grouper );
	//CPPUNIT_TEST( test_3FullLengthFragments );
	//CPPUNIT_TEST( test_3FullLengthFragments_1Truncated );
	//CPPUNIT_TEST(test_runAsScript_bigData);
	//CPPUNIT_TEST(test_runAsScript_smallData_filter_length_40_46);
	//CPPUNIT_TEST(test_runAsScript_loadPath_option_grouper_only);
	//CPPUNIT_TEST(test_runAsScript_loadPath_option);
	//CPPUNIT_TEST(test_runAsScript_loadPath_option_small_data);
 	CPPUNIT_TEST_SUITE_END();

 	SDGString grouperFile;

	public:
		void setUp();
		void tearDown();

	protected:
		void test_3FullLengthFragments( void );
		void test_3FullLengthFragments_1Truncated( void );
		void test_runAsScript_bigData(void);
		void test_runAsScript_smallData_filter_length_40_46(void);
		void test_runAsScript_loadPath_option_grouper_only( void );
		void test_runAsScript_loadPath_option( void );
		void test_runAsScript_loadPath_option_small_data( void );
	private:
		static void write_input_file_for_test_group_load_path(SDGString pathFileName);
		static void write_map_file_for_test_group_load_path(SDGString mapFileName);
		static void write_input_align_file_for_test_group_loadPath_option_small_data(SDGString alignFileName);

};

#endif
