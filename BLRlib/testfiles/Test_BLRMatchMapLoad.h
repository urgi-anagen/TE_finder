#ifndef TEST_BLRMATCHMAPLOAD_H
#define TEST_BLRMATCHMAPLOAD_H

#include <stdlib.h>
#include <map>
#include <list>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "RangePair.h"
#include "BLRMatcherParameter.h"
#include "BLRMatchMap.h"

class Test_BLRMatchMapLoad: public CPPUNIT_NS::TestFixture
{
  CPPUNIT_TEST_SUITE(Test_BLRMatchMapLoad);
  CPPUNIT_TEST(test_readAlign);
  CPPUNIT_TEST(test_load);
  CPPUNIT_TEST(test_readPath);
  CPPUNIT_TEST(test_loadPath);
//  CPPUNIT_TEST(test_loadPath_name2Num_and_num2Name);
//  CPPUNIT_TEST(test_loadPath_rpsList);
//  CPPUNIT_TEST(test_loadPath_with_join_end_of_file);
//  CPPUNIT_TEST(test_loadPath_with_join_begin_of_file);
//  CPPUNIT_TEST(test_loadPath_with_join_middle_of_file);
  CPPUNIT_TEST_SUITE_END();

 public:
    void setUp();
    void tearDown();

 protected:
    void test_readAlign(void);
    void test_readPath(void);
    void test_load(void);
    void test_loadPath(void);
    void test_loadPath_name2Num_and_num2Name(void);
    void test_loadPath_rpsList(void);
    void test_loadPath_with_join_end_of_file(void);
    void test_loadPath_with_join_begin_of_file(void);
    void test_loadPath_with_join_middle_of_file(void);
 
 private:
    bool areTwoRpsEqualsWithoutIndentity(RangePairSet rps1, RangePairSet rps2);

};
#endif
