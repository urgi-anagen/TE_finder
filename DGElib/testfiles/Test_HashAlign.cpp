/*
 * \file Test_HashAlign.cpp
 */

#include <iostream>

#include "Test_HashAlign.h"


CPPUNIT_TEST_SUITE_REGISTRATION(Test_HashAlign);

SDGBioSeqDB Test_HashAlign::initBioSeqDB(void){
	SDGString filename = "./data/seqCluster2.fa";
	SDGBioSeqDB db( filename );
	return db;
}
void Test_HashAlign::test_search_on_loop(void){
	SDGBioSeqDB db = initBioSeqDB();
	/*
	for( SDGBioSeqDB::iterator i=db.begin(); i!=db.end();i++ )
	{
		for( SDGBioSeqDB::iterator j=i+1; j!=db.end();j++)
		{
			HashAlign al;
			al.setSeq(*i,*j);
			al.search();
		}
	}*/
}

