/*
 * \file Test_HashAlign.cpp
 */

#include <iostream>

#include "Test_HashAlignClone.h"


CPPUNIT_TEST_SUITE_REGISTRATION(Test_HashAlignClone);


void Test_HashAlignClone::test_add_push_back_case(void){
	HashAlignClone alc = HashAlignClone();
	
  	std::map<int,std::list<Range> > diag_map;
	int index = 1;
	Range r = Range( 1, 10 );
	alc.addToDiagMap(r, index, diag_map);
	
	Range expRange = Range(1, 10);
	Range obsRange = diag_map[index].back();
	/*	
	std::cout<<" "<<std::endl;	
	std::cout<<"Exp "<<std::endl;
	expRange.view();
	
	std::cout<<" "<<std::endl;	
	std::cout<<"Obs "<<std::endl;
	obsRange.view();
	*/
	CPPUNIT_ASSERT(expRange == obsRange);
	

}

void Test_HashAlignClone::test_add_set_end_case(void)
{
	HashAlignClone alc = HashAlignClone();
	std::map<int,std::list<Range> > diag_map;
	int index = 1;
	
	Range r = Range( 1, 10 );
	alc.addToDiagMap(r, index, diag_map);

	Range r2 = Range(1, 11);
	alc.addToDiagMap(r2, index, diag_map);
	
	Range expRange = Range(1, 11);
	Range obsRange = diag_map[index].back();

	CPPUNIT_ASSERT(expRange == obsRange);
	
	int expSize = 1;
	int obsSize = diag_map[index].size();		

	CPPUNIT_ASSERT(expSize == obsSize);
}

void Test_HashAlignClone::test_add_set_start_case(void)
{
	HashAlignClone alc = HashAlignClone();
	std::map<int,std::list<Range> > diag_map;
	int index = 1;
	
	Range r = Range( 2, 10 );
	alc.addToDiagMap(r, index, diag_map);

	Range r2 = Range(1, 10);
	alc.addToDiagMap(r2, index, diag_map);
	
	Range expRange = Range(1, 10);
	Range obsRange = diag_map[index].back();

	CPPUNIT_ASSERT(expRange == obsRange);
	
	int expSize = 1;
	int obsSize = diag_map[index].size();		

	CPPUNIT_ASSERT(expSize == obsSize);
}
