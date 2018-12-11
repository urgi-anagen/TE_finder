#include "Test_DiagClone.h"
#include "../Range.h"

CPPUNIT_TEST_SUITE_REGISTRATION( Test_DiagClone );

void Test_DiagClone::test_add_push_back_case(void)
{
	DiagClone dc;
	Range r = Range( 1, 10 );
	dc.add(r);

	Range expRange = Range(1, 10);
	Range obsRange = dc.back();

	CPPUNIT_ASSERT(expRange == obsRange);
			
}

void Test_DiagClone::test_add_set_end_case(void)
{
	DiagClone dc;
	Range r = Range( 1, 10 );
	dc.add(r);

	Range r2 = Range(1, 11);
	dc.add(r2);
	
	Range expRange = Range(1, 11);
	Range obsRange = dc.back();

	CPPUNIT_ASSERT(expRange == obsRange);
	
	int expSize = 1;
	int obsSize = dc.size();		

	CPPUNIT_ASSERT(expSize == obsSize);
}

void Test_DiagClone::test_add_set_start_case(void)
{
	DiagClone dc;
	Range r = Range( 2, 10 );
	dc.add(r);

	Range r2 = Range(1, 10);
	dc.add(r2);
	
	Range expRange = Range(1, 10);
	Range obsRange = dc.back();

	CPPUNIT_ASSERT(expRange == obsRange);
	
	int expSize = 1;
	int obsSize = dc.size();		

	CPPUNIT_ASSERT(expSize == obsSize);
}
