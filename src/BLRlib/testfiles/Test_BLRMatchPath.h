#ifndef TEST_BLRMATCHPATH_H
#define TEST_BLRMATCHPATH_H

#include <stdlib.h>
#include <map>
#include <list>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "RangePair.h"
#include "../../matcher/BLRMatcherParameter.h"
#include "BLRMatchMap.h"

class Test_BLRMatchPath: public CPPUNIT_NS::TestFixture
{
  CPPUNIT_TEST_SUITE(Test_BLRMatchPath);
  CPPUNIT_TEST(test_read);
  CPPUNIT_TEST(test_load);
  CPPUNIT_TEST(test_writeBED);
  CPPUNIT_TEST_SUITE_END();

 public:
    void setUp();
    void tearDown();

 protected:
    void test_read(void);
    void test_load(void);
    void test_writeBED(void);

};
#endif
