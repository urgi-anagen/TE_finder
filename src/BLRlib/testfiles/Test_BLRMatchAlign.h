#ifndef TEST_BLRMATCHALIGN_H
#define TEST_BLRMATCHALIGN_H

#include <stdlib.h>
#include <map>
#include <list>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "RangePair.h"
#include "../../matcher/BLRMatcherParameter.h"
#include "BLRMatchMap.h"

class Test_BLRMatchAlign: public CPPUNIT_NS::TestFixture
{
  CPPUNIT_TEST_SUITE(Test_BLRMatchAlign);
  CPPUNIT_TEST(test_read);
  CPPUNIT_TEST(test_load);
  CPPUNIT_TEST(test_set);
  CPPUNIT_TEST_SUITE_END();

    BLRMatcherThreadsParameter createParameter(void) {
        BLRMatcherThreadsParameter para;
        para.setLenFilter(0);
        para.setEvalFilter(10);
        para.setIdFilter(0);
        return para;
    };

 public:
    void setUp();
    void tearDown();

 protected:
    void test_read(void);
    void test_load(void);
    void test_set(void);

};
#endif
