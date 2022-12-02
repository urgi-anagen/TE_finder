//
// Created by Hadi Quesneville on 31/05/2022.
//

#ifndef TE_FINDER_TEST_FASTEXTALIGN_H
#define TE_FINDER_TEST_FASTEXTALIGN_H
#include <cppunit/extensions/HelperMacros.h>

#include "FastExtAlign.h"

class Test_FastExtAlign : public CppUnit::TestFixture{
    CPPUNIT_TEST_SUITE(Test_FastExtAlign);

    CPPUNIT_TEST( test_align );

    CPPUNIT_TEST_SUITE_END();

public:

    Test_FastExtAlign(void) {}

    void setUp()
    {
    }
    void tearDown()
    {
    }

protected:
    void test_align(void);

};


#endif //TE_FINDER_TEST_FASTEXTALIGN_H

