#ifndef TEST_BLRMATCHJOIN_H
#define TEST_BLRMATCHJOIN_H

#include <stdlib.h>
#include <map>
#include <list>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "RangePair.h"
#include "../../matcher/BLRMatcherParameter.h"
#include "BLRMatchMap.h"

class Test_BLRMatchJoin: public CPPUNIT_NS::TestFixture
{
  CPPUNIT_TEST_SUITE(Test_BLRMatchJoin);
  CPPUNIT_TEST(test_join);
  CPPUNIT_TEST(test_nojoin);
  CPPUNIT_TEST(test_clean_conflicts);
  CPPUNIT_TEST(test_merge);
  CPPUNIT_TEST(test_split);
  CPPUNIT_TEST_SUITE_END();

 public:
    void setUp();
    void tearDown();

 protected:
    void test_join(void);
    void test_nojoin(void);
    void test_clean_conflicts(void);
    void test_merge(void);
    void test_split(void);

};
#endif
