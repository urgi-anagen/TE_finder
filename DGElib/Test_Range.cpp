#include "Test_Range.h"

CPPUNIT_TEST_SUITE_REGISTRATION( Test_Range );

void Test_Range::setUp()
{
}

void Test_Range::tearDown()
{
}

void Test_Range::test_getStart( void )
{
	Range i = Range( 1, 10 );
	unsigned long exp = 1;
	unsigned long obs = i.getStart();
	CPPUNIT_ASSERT_EQUAL( obs, exp );
	i = Range( 10, 1 );
	exp = 10;
	obs = i.getStart();
	CPPUNIT_ASSERT_EQUAL( obs, exp );
}

void Test_Range::test_getEnd( void )
{
	Range i = Range( 1, 10 );
	unsigned long exp = 10;
	unsigned long obs = i.getEnd();
	CPPUNIT_ASSERT_EQUAL( obs, exp );
	i = Range( 10, 1 );
	exp = 1;
	obs = i.getEnd();
	CPPUNIT_ASSERT_EQUAL( obs, exp );
}

void Test_Range::test_isPlusStrand( void )
{
	Range i = Range( 1, 10 );
	bool obs = i.isPlusStrand();
	CPPUNIT_ASSERT( obs );
	i = Range( 10, 1 );
	obs = i.isPlusStrand();
	CPPUNIT_ASSERT( ! obs );
}

void Test_Range::test_getLength( void )
{
	Range i = Range( 1, 10 );
	unsigned long exp = 10;
	unsigned long obs = i.getLength();
	CPPUNIT_ASSERT_EQUAL( obs, exp );
	i = Range( 10, 1 );
	obs = i.getLength();
	CPPUNIT_ASSERT_EQUAL( obs, exp );
}

void Test_Range::test_getMin( void )
{
	Range i = Range( 1, 10 );
	unsigned long exp = 1;
	unsigned long obs = i.getMin();
	CPPUNIT_ASSERT_EQUAL( obs, exp );
	i = Range( 10, 1 );
	obs = i.getMin();
	CPPUNIT_ASSERT_EQUAL( obs, exp );
}

void Test_Range::test_getMax( void )
{
	Range i = Range( 1, 10 );
	unsigned long exp = 10;
	unsigned long obs = i.getMax();
	CPPUNIT_ASSERT_EQUAL( obs, exp );
	i = Range( 10, 1 );
	obs = i.getMax();
	CPPUNIT_ASSERT_EQUAL( obs, exp );
}

void Test_Range::test_reverse( void )
{
	Range i = Range( 1, 10 );
	Range exp = Range( 10, 1 );
	i.reverse();
	CPPUNIT_ASSERT( i == exp );
	i = Range( 10, 1 );
	exp = Range( 1, 10 );
	i.reverse();
	CPPUNIT_ASSERT( i == exp );
}

void Test_Range::test_overlap( void )
{
	Range i = Range( 1, 10 );
	Range j = Range( 10, 15 );
	bool obs = i.overlap( j );
	CPPUNIT_ASSERT( obs );
	i = Range( 1, 10 );
	j = Range( 2, 8 );
	obs = i.overlap( j );
	CPPUNIT_ASSERT( obs );
	i = Range( 1, 10 );
	j = Range( 12, 3 );
	obs = i.overlap( j );
	CPPUNIT_ASSERT( obs );
}

void Test_Range::test_isContained( void )
{
	Range i = Range( 2, 8 );
	Range j = Range( 1, 10 );
	bool obs = i.isContained( j );
	CPPUNIT_ASSERT( obs );
	i = Range( 1, 10 );
	j = Range( 1, 10 );
	obs = i.overlap( j );
	CPPUNIT_ASSERT( obs );
}


void Test_Range::test_diff(void)
{
	Range range1 = Range(109951, 110569);
	Range range2 = Range(109985, 110398);
        /*	
	std::cout<<" "<<std::endl;
	std::cout<<"before diff "<<std::endl;
    	std::cout<<"range 1: "<<std::endl;
    	range1.view();
     	std::cout<<"range 2: "<<std::endl;
    	range2.view();
	*/
	Range obsOutRange = range1.diff(range2);
	/*
	std::cout<<"after diff "<<std::endl;
    	std::cout<<"range 1: "<<std::endl;
    	range1.view();
     	std::cout<<"range 2: "<<std::endl;
    	range2.view();
     	std::cout<<"outRange: "<<std::endl;
    	obsOutRange.view();
	*/
	Range expRange1 = Range(109951, 109984);
	Range expOutRange = Range(110398, 110569);
	
	
	CPPUNIT_ASSERT (expRange1 == range1);		
	CPPUNIT_ASSERT (expOutRange == obsOutRange);		
}

void Test_Range::test_diff_range1_start_changed_after_diff(void)
{
	Range range1 = Range(150, 300);
	Range range2 = Range(100, 200);
        /*	
	std::cout<<" "<<std::endl;
	std::cout<<"before diff "<<std::endl;
    	std::cout<<"range 1: "<<std::endl;
    	range1.view();
     	std::cout<<"range 2: "<<std::endl;
    	range2.view();
	*/
	range1.diff(range2);
	/*	
	std::cout<<"after diff "<<std::endl;
    	std::cout<<"range 1: "<<std::endl;
    	range1.view();
     	std::cout<<"range 2: "<<std::endl;
    	range2.view();
     	std::cout<<"outRange: "<<std::endl;
    	obsOutRange.view();
	*/
	Range expRange1 = Range(201, 300);
		
	CPPUNIT_ASSERT (expRange1 == range1);		
}


