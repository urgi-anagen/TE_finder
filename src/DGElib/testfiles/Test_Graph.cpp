//
// Created by Hadi Quesneville on 2019-02-05.
//

#include "Test_Graph.h"

CPPUNIT_TEST_SUITE_REGISTRATION( Test_Graph );

void Test_Graph::setUp()
{
}

void Test_Graph::tearDown()
{
}

void Test_Graph::test_compConnex( void )
{

    Graph<char> gr;

    gr.add_edge('a','b');
    gr.add_edge('b','c');
    gr.add_edge('b','d');

    gr.add_edge('e','f');
    gr.add_edge('f','g');

    gr.add_edge('h','i');

    gr.add_node('j');

    std::ostringstream exp,obs;

    exp<<"*a-b-c-d-\n";
    exp<<"*e-f-g-\n";
    exp<<"*h-i-\n";
    exp<<"*j-\n";

    std::vector<std::vector<char>> vec;
    gr.connexComp(vec);
    for (std::vector<std::vector<char> >::iterator it_vec = vec.begin(); it_vec != vec.end(); it_vec++) {
        obs<<"*";
        for (std::vector<char>::iterator it_vec2 = it_vec->begin(); it_vec2 != it_vec->end(); it_vec2++) {
            obs<<*it_vec2<<"-";
        }
        obs<<std::endl;
    }

    CPPUNIT_ASSERT_EQUAL(exp.str(),obs.str());
}