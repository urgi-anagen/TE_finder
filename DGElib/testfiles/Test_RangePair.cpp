/*
 * \file Test_RangePair.cpp
 */

#include <iostream>

#include "Test_RangePair.h"
#include "../RangePair.h"
#include "../RangeAlign.h"
#include "SDGString.h"
#include "../FileUtils.h"
#include <stdlib.h>

CPPUNIT_TEST_SUITE_REGISTRATION( Test_RangePair );

void Test_RangePair::setUp()
{
}

void Test_RangePair::tearDown()
{
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
