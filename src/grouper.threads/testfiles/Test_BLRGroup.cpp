#include "Test_BLRGroup.h"
#include "Test_BLRGroupUtils.h"
#include "../BLRGroupThreads.h"
#include "../BLRMemIdxBin.h"
#include "BLRGrouperParameter.h"
#include "SDGString.h"
#include "FileUtils.h"

#include "BLRMatchMap.h"
#include "Test_BLRUtils.h"

CPPUNIT_TEST_SUITE_REGISTRATION(Test_BLRGroup);
//------------------------------------------------------------------------------------------------------------
void Test_BLRGroup::test_eraseMember( void )
{
    BLRGrouperParameter para;
    BLRMatchMap mm(para);

    //test
    BLRGroup gr(para,mm,3);

	Member m1(RangeAlignSet(RangeAlign(1,100,200)),1);
	Member m2(RangeAlignSet(RangeAlign(2,100,200)),2);
	GROUPLIST::iterator gr_it=gr.addGroup(m1,m2);
	Member m3(RangeAlignSet(RangeAlign(1,110,210)),3);
	Member m4(RangeAlignSet(RangeAlign(3,100,200)),4);
	gr.addMember(gr_it,m3);
	std::list<unsigned>::iterator i=gr.addMember(gr_it,m4);

	gr.eraseMember(gr_it,i);

	//expectation
	BLRGroup gr_exp(para,mm,3);
	Member m1e=m1,m2e=m2,m3e=m3;
	GROUPLIST::iterator gr_exp_it=gr_exp.addGroup(m1e,m2e);
	gr_exp.addMember(gr_exp_it,m3e);

	//compare test to expectation
	std::ostringstream ostr_obs;
	gr.show_group(ostr_obs);
	std::ostringstream ostr_exp;
	gr_exp.show_group(ostr_exp);

	CPPUNIT_ASSERT_EQUAL(ostr_exp.str(), ostr_obs.str());
}

/*void Test_BLRGroup::test_mergeWithExtGroup(void)
{
    para_ptr=new BLRGrouperParameter();
    BLRMatchMap mm1(para_ptr),mm2(para_ptr);


    //Test
    BLRGroup ext_gr(para_ptr,&mm1,3);
    BLRGroup gr(para_ptr,&mm2,3);


    //group init
	Member m11(RangeAlignSet(RangeAlign(1,100,200)),1);
	Member m12(RangeAlignSet(RangeAlign(2,100,200)),2);
	GROUPLIST::iterator gr_it=gr.addGroup(m11,m12);
	Member m13(RangeAlignSet(RangeAlign(1,110,210)),3);
	Member m14(RangeAlignSet(RangeAlign(3,100,200)),4);
	gr.addMember(gr_it,m13);
	gr.addMember(gr_it,m14);

	gr.view_a_group(gr_it);

	//ext group init
	Member m21(RangeAlignSet(RangeAlign(1,100,200)),5);
	Member m22(RangeAlignSet(RangeAlign(2,1000,1200)),6);
	GROUPLIST::iterator ext_gr_it=ext_gr.addGroup(m21,m22);
	Member m23(RangeAlignSet(RangeAlign(1,2000,2100)),7);
	Member m24(RangeAlignSet(RangeAlign(3,1000,1200)),8);
	ext_gr.addMember(ext_gr_it,m23);
	ext_gr.addMember(ext_gr_it,m24);

	ext_gr.view_a_group(ext_gr_it);


	gr.mergeWithExtGroup(gr_it, ext_gr_it, ext_gr);


	//Expectation. Note insertion order is important
    BLRGroup gr_exp(para_ptr,&mm2,3);
	GROUPLIST::iterator gr_exp_it=gr_exp.addGroup(m11,m12);
	gr_exp.addMember(gr_exp_it,m13);
	gr_exp.addMember(gr_exp_it,m14);
	gr_exp.addMember(gr_exp_it,m24);
	gr_exp.addMember(gr_exp_it,m23);
	gr_exp.addMember(gr_exp_it,m22);
	gr_exp.addMember(gr_exp_it,m21);

	//compare test to expectation
	std::ostringstream ostr_obs;
	gr.show_group(ostr_obs);
	std::ostringstream ostr_exp;
	gr_exp.show_group(ostr_exp);

	CPPUNIT_ASSERT_EQUAL(ostr_exp.str(), ostr_obs.str());
}*/


