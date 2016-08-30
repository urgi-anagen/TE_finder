/*
 * \file Test_RangePair.cpp
 */

#include <iostream>

#include "Test_RangePair.h"
#include "RangePair.h"
#include "RangeAlign.h"
#include "SDGString.h"
#include "FileUtils.h"
#include <stdlib.h>

CPPUNIT_TEST_SUITE_REGISTRATION( Test_RangePair );

void Test_RangePair::setUp()
{
}

void Test_RangePair::tearDown()
{
}

void Test_RangePair::test_set_string( void )
{
	RangePair rangePairObs = RangePair( "seq1\t1\t100\tseq2\t1\t100\t1e-187\t213\t97.2" );

	RangePair rangePairExp = RangePair();
	RangeAlign rangeQuery = RangeAlign( "seq1", 0, 1, 100 );
	RangeAlign rangeSubject = RangeAlign( "seq2", 0, 1, 100 );
	rangePairExp.setRangeQ( rangeQuery );
	rangePairExp.setRangeS( rangeSubject );
	rangePairExp.setE_value( 1e-187 );
	rangePairExp.setScore( 213 );
	rangePairExp.setIdentity(97.2);

	bool exp = true;
	bool obs = (rangePairExp == rangePairObs);

	CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangePair::test_writetxt( void )
{
	SDGString expFile = "dummyExpFile";
	std::ofstream expFileStream;
	FileUtils::openFile( expFile, expFileStream );
	expFileStream<<"seq1\t1\t100\tseq2\t1\t100\t1e-187\t213\t97.2"<<std::endl;
	expFileStream.close();

	SDGString obsFile = "dummyObsFile";
	std::ofstream obsFileStream;
	FileUtils::openFile( obsFile, obsFileStream );
	RangePair rp = RangePair( "seq1\t1\t100\tseq2\t1\t100\t1e-187\t213\t97.2" );
	rp.writetxt( obsFileStream );
	obsFileStream.close();

	bool exp = true;
	bool obs = FileUtils::areTwoFilesIdentical( expFile, obsFile );
	CPPUNIT_ASSERT_EQUAL( exp, obs );

	FileUtils::openFile( obsFile, obsFileStream );
	rp.set( "seq2\t1\t100\tseq1\t1\t100\t1e-187\t213\t97.2" );
	rp.writetxt( obsFileStream );
	obsFileStream.close();

	exp = false;
	obs = FileUtils::areTwoFilesIdentical( expFile, obsFile );
	CPPUNIT_ASSERT_EQUAL( exp, obs );

	FileUtils::removeFile( expFile );
	FileUtils::removeFile( obsFile );
}

void Test_RangePair::test_diffQ(void)
{
	RangePair rp1 = RangePair("1\t109951\t110569\t2\t119\t689\t9.4e-19\t619\t77.768");
	RangePair rp2 = RangePair("1\t109985\t110398\t1\t2429\t2033\t0\t414\t77.8604");

	RangePair obsRpOut = rp1.diffQ(rp2);
	RangePair obsRp1 = rp1;
     
	RangePair expRp1 = RangePair("1\t109951\t109984\t2\t119\t149\t9.4e-19\t34\t77.768");
	RangePair expRpOut = RangePair("1\t110398\t110569\t2\t531\t689\t9.4e-19\t172\t77.768");
	/*	
    	std::cout<<" "<<std::endl;
    	std::cout<<"obs rp 1: "<<std::endl;
    	obsRp1.view();
		
	std::cout<<" "<<std::endl;
    	std::cout<<"exp rp 1: "<<std::endl;
    	expRp1.view();
	*/
 	CPPUNIT_ASSERT (expRp1 == obsRp1);
	CPPUNIT_ASSERT (expRpOut == obsRpOut);	
}

void Test_RangePair::test_diffQ_rp1_start_modified(void)
{
	RangePair rp1 = RangePair("1\t150\t300\t2\t119\t689\t9.4e-19\t619\t77.768");
	RangePair rp2 = RangePair("1\t100\t200\t1\t2429\t2033\t0\t414\t77.8604");

	RangePair obsRpOut = rp1.diffQ(rp2);
	RangePair obsRp1 = rp1;
     
	RangePair expRp1 = RangePair("1\t201\t300\t2\t312\t689\t9.4e-19\t410\t77.768");
	/*	
    	std::cout<<" "<<std::endl;
    	std::cout<<"obs rp 1: "<<std::endl;
    	obsRp1.view();
		
	std::cout<<" "<<std::endl;
    	std::cout<<"exp rp 1: "<<std::endl;
    	expRp1.view();
	*/
 	CPPUNIT_ASSERT (expRp1 == obsRp1);
}
/*
void Test_RangePair::test_diffQ_rp1_start_modified_reverse_strand(void)
{
	RangePair rp1 = RangePair("1\t150\t300\t2\t119\t689\t9.4e-19\t619\t77.768");
	RangePair rp2 = RangePair("1\t100\t200\t1\t2429\t2033\t0\t414\t77.8604");

	RangePair obsRpOut = rp1.diffQ(rp2);
	RangePair obsRp1 = rp1;
     
	RangePair expRp1 = RangePair("1\t201\t300\t2\t312\t689\t9.4e-19\t410\t77.768");
		

    	std::cout<<" "<<std::endl;
    	std::cout<<"obs rp 1: "<<std::endl;
    	obsRp1.view();
		
	std::cout<<" "<<std::endl;
    	std::cout<<"exp rp 1: "<<std::endl;
    	expRp1.view();
	
 	CPPUNIT_ASSERT (expRp1 == obsRp1);
}
*/

void Test_RangePair::test_reComputeSubjectCoords(void)
{
	RangePair rp1 = RangePair("1\t109951\t110569\t2\t119\t689\t9.4e-19\t619\t77.768");
	RangePair rp2 = RangePair("1\t109985\t110398\t1\t2429\t2033\t0\t414\t77.8604");

	unsigned queryStart = rp1.getRangeQ().getStart(); 
	unsigned queryEnd = rp1.getRangeQ().getEnd();

	RangeAlign rpOut = rp1.getRangeQ().diff(rp2.getRangeQ());

	RangePair newRp = rp1;
	newRp.setRangeQ(rpOut);

	rp1.reComputeSubjectCoords(newRp, queryStart, queryEnd);
	
	unsigned long expSubjectStartOnRp1 = 119;
	unsigned long expSubjectEndOnRp1 = 149;
	
	unsigned long expSubjectStartOnNewRp = 531;
	unsigned long expSubjectEndOnNewRp = 689;

	unsigned long obsSubjectStartOnRp1 = rp1.getRangeS().getStart();
	unsigned long obsSubjectEndOnRp1 = rp1.getRangeS().getEnd();

	unsigned long obsSubjectStartOnNewRp = newRp.getRangeS().getStart();
	unsigned long obsSubjectEndOnNewRp = newRp.getRangeS().getEnd();

	/*		 
	std::cout<<" "<<std::endl;
    	std::cout<<"rp1 "<<std::endl;
    	rp1.view();
	
	std::cout<<" "<<std::endl;
    	std::cout<<"newRp "<<std::endl;
    	newRp.view();
	*/
	CPPUNIT_ASSERT (expSubjectStartOnRp1 == obsSubjectStartOnRp1);	
	CPPUNIT_ASSERT (expSubjectEndOnRp1 == obsSubjectEndOnRp1);
	CPPUNIT_ASSERT (expSubjectStartOnNewRp == obsSubjectStartOnNewRp);
	CPPUNIT_ASSERT (expSubjectEndOnNewRp == obsSubjectEndOnNewRp);	
}

void Test_RangePair::test_reComputeSubjectCoords_reverse_strand(void)
{
	RangePair rp1 = RangePair("1\t109951\t110569\t2\t689\t119\t9.4e-19\t619\t77.768");
	RangePair rp2 = RangePair("1\t110398\t109985\t1\t2429\t2033\t0\t414\t77.8604");

	unsigned queryStart = rp1.getRangeQ().getStart(); 
	unsigned queryEnd = rp1.getRangeQ().getEnd();

	RangeAlign rpOut = rp1.getRangeQ().diff(rp2.getRangeQ());

	RangePair newRp = rp1;
	newRp.setRangeQ(rpOut);

	rp1.reComputeSubjectCoords(newRp, queryStart, queryEnd);
	
	unsigned long expSubjectStartOnRp1 = 689;
	unsigned long expSubjectEndOnRp1 = 659;
	
	unsigned long expSubjectStartOnNewRp = 277;
	unsigned long expSubjectEndOnNewRp = 119;

	unsigned long obsSubjectStartOnRp1 = rp1.getRangeS().getStart();
	unsigned long obsSubjectEndOnRp1 = rp1.getRangeS().getEnd();

	unsigned long obsSubjectStartOnNewRp = newRp.getRangeS().getStart();
	unsigned long obsSubjectEndOnNewRp = newRp.getRangeS().getEnd();

	/*			 
	std::cout<<" "<<std::endl;
    	std::cout<<"rp1 "<<std::endl;
    	rp1.view();
	
	std::cout<<" "<<std::endl;
    	std::cout<<"newRp "<<std::endl;
    	newRp.view();
	*/
	CPPUNIT_ASSERT (expSubjectStartOnRp1 == obsSubjectStartOnRp1);	
	CPPUNIT_ASSERT (expSubjectEndOnRp1 == obsSubjectEndOnRp1);
	CPPUNIT_ASSERT (expSubjectStartOnNewRp == obsSubjectStartOnNewRp);
	CPPUNIT_ASSERT (expSubjectEndOnNewRp == obsSubjectEndOnNewRp);	
}
void Test_RangePair::test_getRangeQ_affectation(void)
{
	RangePair rp1 = RangePair("1\t109951\t110569\t2\t119\t689\t9.4e-19\t619\t77.768");
	RangeAlign rangeQuery = RangeAlign( "seq1", 0, 1, 100 );
	/*	 
	std::cout<<" "<<std::endl;
    	std::cout<<"before rp1.getRangeQ() before affectation "<<std::endl;
    	rp1.getRangeQ().view();
	*/
	rp1.getRangeQ() = rangeQuery;
	RangeAlign obsRangeQ = rp1.getRangeQ();
	/*  
	std::cout<<" "<<std::endl;
    	std::cout<<"before rp1.getRangeQ() after affectation "<<std::endl;
    	obsRangeQ.view();
	*/
	RangeAlign expRangeQ = RangeAlign("seq1", 0, 1, 100);	
	CPPUNIT_ASSERT (expRangeQ == obsRangeQ); 
}

void Test_RangePair::test_isPlusStrand_reverse_strand(void)
{

	RangePair rp1 = RangePair("1\t109951\t110569\t2\t689\t119\t9.4e-19\t619\t77.768");
	bool obsIsPlusStrand = rp1.isPlusStrand();
	bool expIsPluqstrand = false;

	CPPUNIT_ASSERT_EQUAL(obsIsPlusStrand, expIsPluqstrand);
}
/*
void Test_RangePair::test_mergeQ(void)
{
	RangePair rp1 = test_createRangePairSet1ForTest_mergeQ();	
	RangePair rp2 = test_createRangePairSet2ForTest_mergeQ();	
	
	std::cout<<"rps1 test_mergeQ"<<std::endl;
	rp1.view();

	std::cout<<"rps2 test_mergeQ"<<std::endl;
	rp2.view();


}

RangePair Test_RangePair::test_createRangePairSet1ForTest_mergeQ(void)
{
  std::list<RangePair> rpList2;
  SDGString line2 = " \t105085\t110398\t \t3202\t2033\t0\t683\t80.6808";
  RangePairSet rangePairSet2 = RangePairSet(line2);
  rangePairSet2.getRangeQ().setNumChr(1);
  rangePairSet2.getRangeQ().setNameSeq("");
  rangePairSet2.getRangeS().setNumChr(1);
  rangePairSet2.getRangeS().setNameSeq("");
  rangePairSet2.setLength(683);
  // match part   
  SDGString line21 = " \t109985\t110398\t \t2429\t2033\t0\t414\t77.8604";
  RangePair rangePair21 = RangePair(line21);
  rangePair21.getRangeQ().setNumChr(1);
  rangePair21.getRangeQ().setNameSeq("");
  rangePair21.getRangeS().setNumChr(1);
  rangePair21.getRangeS().setNameSeq("");
  rangePair21.setLength(414);
  rpList2.push_back(rangePair21);
  // match part   
  SDGString line22 = " \t105085\t105353\t \t3202\t2921\t0\t269\t84.8361";
  RangePair rangePair22 = RangePair(line22);
  rangePair22.getRangeQ().setNumChr(1);
  rangePair22.getRangeQ().setNameSeq("");
  rangePair22.getRangeS().setNumChr(1);
  rangePair22.getRangeS().setNameSeq("");
  rangePair22.setLength(269);
  // set match part to match
  rpList2.push_back(rangePair22);

  rangePairSet2.setPathDirectly(rpList2);

  return rangePairSet2;
}

RangePair Test_RangePair::test_createRangePairSet2ForTest_mergeQ(void)
{
     // match 
  SDGString line1 = " \t109951\t110569\t \t119\t689\t9.4e-19\t619\t77.768";
  RangePairSet rangePairSet1 = RangePairSet(line1);
  rangePairSet1.getRangeQ().setNumChr(1);
  rangePairSet1.getRangeQ().setNameSeq("");
  rangePairSet1.getRangeS().setNumChr(2);
  rangePairSet1.getRangeS().setNameSeq("");
  rangePairSet1.setLength(619);
  // match part
  SDGString line11 = " \t109951\t110569\t \t119\t689\t9.4e-19\t619\t77.768";
  RangePair rangePair11 = RangePair(line11);
  rangePair11.getRangeQ().setNumChr(1);
  rangePair11.getRangeQ().setNameSeq("");
  rangePair11.getRangeS().setNumChr(2);
  rangePair11.getRangeS().setNameSeq("");
  rangePair11.setLength(619);

  // set match part to match
  std::list<RangePair> rpList1;
  rpList1.push_back(rangePair11);

  rangePairSet1.setPathDirectly(rpList1);
  return rangePairSet1;


}
*/

