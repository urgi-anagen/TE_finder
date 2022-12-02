#include <SDGMemBioSeq.h>

#include "Test_FragAlign.h"

CPPUNIT_TEST_SUITE_REGISTRATION(Test_FragAlign);
//------------------------------------------------------------------------------------------------------------
void Test_FragAlign::test_join( void )
{

    std::list< RangePair > frag;

    frag.push_back(RangePair(1, 10, 50,1,110,150,40,0.0, 0.0));
    frag.push_back(RangePair(1, 60, 150,1,160,200,40,0.0, 0.0));
    frag.push_back(RangePair(1, 10, 40,1,120,140,40,0.0, 0.0));
    frag.push_back(RangePair(1, 10, 50,2,110,150,40,0.0, 0.0));
    frag.push_back(RangePair(1, 60, 150,2,160,200,40,0.0, 0.0));
    frag.push_back(RangePair(1, 10, 40,2,120,140,40,0.0, 0.0));

    double gap_pen=0.5;
    double dist_pen=0.5;
    FragAlign fragAlign(dist_pen, 0, gap_pen,0);
    std::list< RangePairSet > jfrag=fragAlign.join(frag);


    jfrag.sort(RangePairSet::less);
    std::ostringstream ostr_obs;
    for(std::list< RangePairSet >::iterator rp_it=jfrag.begin(); rp_it!=jfrag.end();rp_it++)
    {
        rp_it->write_raw(ostr_obs);
    }

    std::list< RangePair > frag_exp;

    frag_exp.push_back(RangePair(1, 10, 150,1,110,200,75,0.0, 0.0));
    frag_exp.push_back(RangePair(1, 10, 40,1,120,140,40,0.0, 0.0));
    frag_exp.push_back(RangePair(1, 10, 150,2,110,200,75,0.0, 0.0));
    frag_exp.push_back(RangePair(1, 10, 40,2,120,140,40,0.0, 0.0));

    frag_exp.sort(RangePair::less);
    std::ostringstream ostr_exp;
    for(std::list< RangePair >::iterator rp_it=frag_exp.begin(); rp_it!=frag_exp.end();rp_it++)
    {
        rp_it->write_raw(ostr_exp);
    }

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());

}