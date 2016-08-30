/*
 * Test_BLRGroup.h
 *
 *  Created on: 30 déc. 2009
 *      Author: hquesnev
 */

#ifndef TEST_BLRGROUP_H_
#define TEST_BLRGROUP_H_
#include <cppunit/extensions/HelperMacros.h>
#include "Range.h"
#include "RangeAlign.h"
#include "RangeAlignSet.h"
#include "../BLRGroup.h"
#include <list>
#include <vector>

class Test_BLRGroup : public CppUnit::TestFixture {

	BLRGroup* p_blrgroup;
	BLRGrouperParameter* p_para;

	CPPUNIT_TEST_SUITE(Test_BLRGroup);
	//CPPUNIT_TEST(test_testAlignMemberOverlap);
	//CPPUNIT_TEST( test_group );
	//CPPUNIT_TEST( test_group_load_path );
	//CPPUNIT_TEST( test_compare_join_for_debug );
	//CPPUNIT_TEST( test_group_compare_intermediate_rpsList_between_2_loads );
	
	CPPUNIT_TEST_SUITE_END();

public:

	Test_BLRGroup(void) : p_blrgroup(0), p_para(0){}

	void setUp()
	{
	    p_para=new BLRGrouperParameter();

	    int numarg=5;
	    std::vector<SDGString> tabarg(numarg);
	    tabarg[0]="./test";
	    tabarg[1]="-m Harm_BAC.align";
	    tabarg[2]="-q Harm_BAC.fa";
	    tabarg[3]="-C 0.95";
	    tabarg[4]="-v 1";

//	    p_para->parseOptArg(numarg,tabarg); // ne trouve pas Harm_BAC.fa dans le répertoire courant
	    // où est lancé le programme. je ne me l'explique pas encore, mais cela ne l'empêche pas de faire
	    // les tests déjà écrits.

	    p_blrgroup=new BLRGroup(p_para);
	}
	void tearDown()
	{
		delete p_blrgroup;
	}

protected:

	void test_testAlignMemberOverlap(void)
	{
		std::list<Range> rs1;
		rs1.insert(rs1.begin(),Range(100,150));
		rs1.insert(rs1.begin(),Range(160,200));
		RangeAlign ra1(1,100,200);
		RangeAlignSet align(ra1,rs1);

		std::list<Range> rs2;
		rs2.insert(rs2.begin(),Range(110,150));
		rs2.insert(rs2.begin(),Range(160,190));
		RangeAlign ra2(1,110,190);
		Member member1(ra2);

		//bool obs=p_blrgroup->testAlignMemberOverlap(align,member1);
		//bool exp=true;
		//CPPUNIT_ASSERT_EQUAL(obs,exp);

		std::list<Range> rs3;
		rs3.insert(rs3.begin(),Range(10,50));
		rs3.insert(rs3.begin(),Range(60,90));
		RangeAlign ra3(1,10,90);
		Member member2(ra3);

		//obs=p_blrgroup->testAlignMemberOverlap(align,member2);
		//exp=false;
		//CPPUNIT_ASSERT_EQUAL(obs,exp);

	}
	
	void test_group(void);
	void test_group_load_path(void);
	void test_compare_join_for_debug(void);
	void test_group_compare_intermediate_rpsList_between_2_loads(void);

};
CPPUNIT_TEST_SUITE_REGISTRATION(Test_BLRGroup);

#endif /* TEST_BLRGROUP_H_ */
