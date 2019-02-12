//
// Created by Hadi Quesneville on 2019-02-05.
//

#ifndef TE_FINDER_TEST_GRAPH_H
#define TE_FINDER_TEST_GRAPH_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../Graph.h"

class Test_Graph : public CPPUNIT_NS::TestFixture {

CPPUNIT_TEST_SUITE(Test_Graph);

CPPUNIT_TEST(test_compConnex);

CPPUNIT_TEST_SUITE_END();

public:
    void setUp();

    void tearDown();

protected:
    void test_compConnex(void);

};


#endif //TE_FINDER_TEST_GRAPH_H
