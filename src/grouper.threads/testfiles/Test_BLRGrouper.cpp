//
// Created by Hadi Quesneville on 2019-02-14.
//

#include "Test_BLRGrouper.h"

CPPUNIT_TEST_SUITE_REGISTRATION(Test_BLRGrouper);
void Test_BLRGrouper::test_mergeGroupsLists(void)
{
    BLRGrouperParameter para;
    BLRMatchMap mm1(para),mm2(para),mm3(para);
    BLRGrouper grpr(para);

    //Test
    BLRGroup ext_gr(para, mm1, 3);
    BLRGroup gr(para, mm2, 3);

    //group init
    Member m11(RangeAlignSet(RangeAlign(1,100,200)),1);
    Member m12(RangeAlignSet(RangeAlign(2,100,200)),2);
    GROUPLIST::iterator gr_it=gr.addGroup(m11,m12);
    Member m13(RangeAlignSet(RangeAlign(1,110,210)),3);
    Member m14(RangeAlignSet(RangeAlign(3,100,200)),4);
    gr.addMember(gr_it,m13);
    gr.addMember(gr_it,m14);

    //ext group init
    Member m21(RangeAlignSet(RangeAlign(1,100,200)),5);
    Member m22(RangeAlignSet(RangeAlign(2,1000,1200)),6);
    GROUPLIST::iterator ext_gr_it=ext_gr.addGroup(m21,m22);
    Member m23(RangeAlignSet(RangeAlign(1,2000,2100)),7);
    Member m24(RangeAlignSet(RangeAlign(3,1000,1200)),8);
    ext_gr.addMember(ext_gr_it,m23);
    ext_gr.addMember(ext_gr_it,m24);

    grpr.mergeGroupsLists(&gr, &ext_gr,0);

    //Expectation. Note insertion order is important
    BLRGroup gr_exp(para, mm3, 3);
    m11.merge(m21);
    m11.idlist.splice( m11.idlist.end(), m21.idlist );

    GROUPLIST::iterator gr_exp_it=gr_exp.addGroup(m11,m12);
    gr_exp.addMember(gr_exp_it,m13);
    gr_exp.addMember(gr_exp_it,m14);
    gr_exp.addMember(gr_exp_it,m24);
    gr_exp.addMember(gr_exp_it,m23);
    gr_exp.addMember(gr_exp_it,m22);

    //compare test to expectation
    std::ostringstream ostr_obs;
    gr.show_group(ostr_obs);
    std::ostringstream ostr_exp;
    gr_exp.show_group(ostr_exp);

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(), ostr_obs.str());
}