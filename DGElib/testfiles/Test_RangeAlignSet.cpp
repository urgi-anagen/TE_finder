#include "Test_RangeAlignSet.h"

CPPUNIT_TEST_SUITE_REGISTRATION( Test_RangeAlignSet );

void Test_RangeAlignSet::setUp()
{
}

void Test_RangeAlignSet::tearDown()
{
}

void Test_RangeAlignSet::test_reset( void )
{
	RangeAlign ra_1 = RangeAlign();
	std::list<Range> lRanges_1;
	lRanges_1.push_back(Range(1, 100));
	lRanges_1.push_back(Range(101, 200));
	RangeAlignSet ras_1 = RangeAlignSet(ra_1, lRanges_1);
	ras_1.sortUsingStrand();
	ras_1.setIncluded(true);

	RangeAlignSet ras_2 = RangeAlignSet();

	ras_1.reset();

	bool exp = true;
	bool obs = (ras_1 == ras_2);
	CPPUNIT_ASSERT_EQUAL(exp, obs);
}

void Test_RangeAlignSet::test_operator_equal_true( void )
{
	RangeAlign ra_1 = RangeAlign();
	std::list<Range> lRanges_1;
	lRanges_1.push_back( Range( 1, 100 ) );
	lRanges_1.push_back( Range( 101, 200 ) );
	RangeAlignSet ras_1 = RangeAlignSet( ra_1, lRanges_1 );
	ras_1.sortUsingStrand();
	ras_1.setIncluded( true );

	RangeAlign ra_2 = RangeAlign();
	std::list<Range> lRanges_2;
	lRanges_2.push_back( Range( 1, 100 ) );
	lRanges_2.push_back( Range( 101, 200 ) );
	RangeAlignSet ras_2 = RangeAlignSet( ra_2, lRanges_2 );
	ras_2.sortUsingStrand();
	ras_2.setIncluded( true );

	bool exp = true;
	bool obs = ( ras_1 == ras_2 );
	CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_operator_equal_false_coordinates( void )
{
	RangeAlign ra_1 = RangeAlign();
	std::list<Range> lRanges_1;
	lRanges_1.push_back( Range( 1, 100 ) );
	lRanges_1.push_back( Range( 101, 200 ) );
	RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );
  ras_1.sortUsingStrand();
  ras_1.setIncluded( true );

  RangeAlign ra_2 = RangeAlign();
  std::list<Range> lRanges_2;
  lRanges_2.push_back( Range( 154, 321 ) );
  lRanges_2.push_back( Range( 584, 799 ) );
  RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );
  ras_2.sortUsingStrand();
  ras_2.setIncluded( true );

  bool exp = false;
  bool obs = ( ras_1 == ras_2 );
  CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_operator_equal_false_included( void )
{
	RangeAlign ra_1 = RangeAlign();
	std::list<Range> lRanges_1;
	lRanges_1.push_back( Range( 1, 100 ) );
	lRanges_1.push_back( Range( 101, 200 ) );
	RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );
  ras_1.sortUsingStrand();
  ras_1.setIncluded( true );

  RangeAlign ra_2 = RangeAlign();
  std::list<Range> lRanges_2;
  lRanges_2.push_back( Range( 1, 100 ) );
  lRanges_2.push_back( Range( 101, 200 ) );
  RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );
  ras_2.sortUsingStrand();
  ras_2.setIncluded( false );

  bool exp = false;
  bool obs = ( ras_1 == ras_2 );
  CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_getStart_plusStrand( void )
{
	RangeAlign ra = RangeAlign();
	std::list<Range> lRanges;
	lRanges.push_back( Range( 1, 100 ) );
  lRanges.push_back( Range( 101, 200 ) );
	RangeAlignSet ras = RangeAlignSet( &ra, lRanges );

	unsigned long exp = 1;
	unsigned long obs = ras.getStart();
	CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_getStart_minusStrand( void )
{
	RangeAlign ra = RangeAlign();
	std::list<Range> lRanges;
	lRanges.push_back( Range( 100, 1 ) );
  lRanges.push_back( Range( 200, 101 ) );
	RangeAlignSet ras = RangeAlignSet( &ra, lRanges );

	unsigned long exp = 200;
	unsigned long obs = ras.getStart();
	CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_getEnd_plusStrand( void )
{
	RangeAlign ra = RangeAlign();
	std::list<Range> lRanges;
	lRanges.push_back( Range( 1, 100 ) );
  lRanges.push_back( Range( 101, 200 ) );
	RangeAlignSet ras = RangeAlignSet( &ra, lRanges );

	unsigned long exp = 200;
	unsigned long obs = ras.getEnd();
	CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_getEnd_minusStrand( void )
{
	RangeAlign ra = RangeAlign();
	std::list<Range> lRanges;
	lRanges.push_back( Range( 100, 1 ) );
  lRanges.push_back( Range( 200, 101 ) );
	RangeAlignSet ras = RangeAlignSet( &ra, lRanges );

	unsigned long exp = 1;
	unsigned long obs = ras.getEnd();
	CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_isPlusStrand_true( void )
{
	RangeAlign ra = RangeAlign();
	std::list<Range> lRanges;
	lRanges.push_back( Range( 1, 100 ) );
  lRanges.push_back( Range( 101, 200 ) );
	RangeAlignSet ras = RangeAlignSet( &ra, lRanges );

	bool exp = true;
	bool obs = ras.isPlusStrand();
	CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_isPlusStrand_false( void )
{
	RangeAlign ra = RangeAlign();
	std::list<Range> lRanges;
	lRanges.push_back( Range( 100, 1 ) );
	lRanges.push_back( Range( 200, 101 ) );
	RangeAlignSet ras = RangeAlignSet( &ra, lRanges );

	bool exp = false;
	bool obs = ras.isPlusStrand();
	CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_sortUsingStrand_plusStrand( void )
{
	RangeAlign ra = RangeAlign();
	std::list<Range> lRanges;
  lRanges.push_back( Range( 101, 200 ) );
	lRanges.push_back( Range( 1, 100 ) );
	RangeAlignSet ras = RangeAlignSet( &ra, lRanges );

  RangeAlign ra_exp = RangeAlign();
  std::list<Range> lRanges_exp;
  lRanges_exp.push_back( Range( 1, 100 ) );
  lRanges_exp.push_back( Range( 101, 200 ) );
  RangeAlignSet ras_exp = RangeAlignSet( &ra_exp, lRanges_exp );

  ras.sortUsingStrand();

  bool exp = true;
  bool obs = ( ras == ras_exp );
  CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_sortUsingStrand_minusStrand( void )
{
	RangeAlign ra = RangeAlign();
	std::list<Range> lRanges;
  lRanges.push_back( Range( 100, 1 ) );
	lRanges.push_back( Range( 200, 101 ) );
	RangeAlignSet ras = RangeAlignSet( &ra, lRanges );

  RangeAlign ra_exp = RangeAlign();
  std::list<Range> lRanges_exp;
  lRanges_exp.push_back( Range( 200, 101 ) );
  lRanges_exp.push_back( Range( 100, 1 ) );
  RangeAlignSet ras_exp = RangeAlignSet( &ra_exp, lRanges_exp );

  ras.sortUsingStrand();

  bool exp = true;
  bool obs = ( ras == ras_exp );
  CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_getIncluded( void )
{
	RangeAlignSet i = RangeAlignSet();
	bool exp = false;
	bool obs = i.getIncluded();
	CPPUNIT_ASSERT_EQUAL( obs, exp );
	exp = true;
	i.setIncluded( true );
	obs = i.getIncluded();
	CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_getLengthSet( void )
{
	RangeAlign ra = RangeAlign();
	std::list<Range> lRanges;
	lRanges.push_back( Range( 1, 100 ) );
	lRanges.push_back( Range( 151, 200 ) );
	RangeAlignSet i = RangeAlignSet( &ra, lRanges );
	unsigned exp = 150;
	unsigned obs = i.getLengthSet();
	CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_reverse( void )
{
	RangeAlign ra_1 = RangeAlign();
	std::list<Range> lRanges_1;
	lRanges_1.push_back( Range( 200, 101 ) );
	lRanges_1.push_back( Range( 100, 1 ) );
	RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );

	RangeAlign ra_exp = RangeAlign();
	std::list<Range> lRanges_exp;
	lRanges_exp.push_back( Range( 1, 100 ) );
	lRanges_exp.push_back( Range( 101, 200 ) );
	RangeAlignSet ras_exp = RangeAlignSet( &ra_exp, lRanges_exp );

	ras_1.reverse();

	CPPUNIT_ASSERT_EQUAL( ras_exp, ras_1 );
}

void Test_RangeAlignSet::test_overlap_length( void )
{
	RangeAlign ra_1 = RangeAlign();
	std::list<Range> lRanges_1;
	lRanges_1.push_back( Range( 1, 100 ) );
	lRanges_1.push_back( Range( 101, 200 ) );
	RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );

	RangeAlign ra_2 = RangeAlign();
	std::list<Range> lRanges_2;
	lRanges_2.push_back( Range( 151, 200 ) );
	lRanges_2.push_back( Range( 301, 400 ) );
	RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );

	unsigned exp = 50;
	unsigned obs = ras_1.overlap_length( ras_2 );
	CPPUNIT_ASSERT_EQUAL( exp, obs );

	lRanges_2.clear();
	lRanges_2.push_back( Range( 200, 151 ) );
	lRanges_2.push_back( Range( 400, 301 ) );
	ras_2 = RangeAlignSet( &ra_2, lRanges_2 );

	exp = 50;
	obs = ras_1.overlap_length( ras_2 );
	CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_isStrictlyContained_sameSize_true( void )
{
	RangeAlign ra_1 = RangeAlign();
	std::list<Range> lRanges_1;
	lRanges_1.push_back( Range( 21, 90 ) );
	lRanges_1.push_back( Range( 155, 190 ) );
	RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );

	RangeAlign ra_2 = RangeAlign();
	std::list<Range> lRanges_2;
	lRanges_2.push_back( Range( 1, 100 ) );
	lRanges_2.push_back( Range( 151, 200 ) );
	RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );

	bool exp = true;
	bool obs = ras_1.isStrictlyContained( ras_2 );
	CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_isStrictlyContained_sameSize_false( void )
{
	RangeAlign ra_1 = RangeAlign();
	std::list<Range> lRanges_1;
	lRanges_1.push_back( Range( 21, 90 ) );
	lRanges_1.push_back( Range( 151, 200 ) );
	RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );

	RangeAlign ra_2 = RangeAlign();
	std::list<Range> lRanges_2;
	lRanges_2.push_back( Range( 1, 100 ) );
	lRanges_2.push_back( Range( 155, 190 ) );
	RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );

	bool exp = false;
	bool obs = ras_1.isStrictlyContained( ras_2 );
	CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_isStrictlyContained_differentSizes_true_2( void )
{
	RangeAlign ra_1 = RangeAlign();
	std::list<Range> lRanges_1;
	lRanges_1.push_back( Range( 21, 90 ) );
	lRanges_1.push_back( Range( 155, 190 ) );
	RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );

	RangeAlign ra_2 = RangeAlign();
	std::list<Range> lRanges_2;
	lRanges_2.push_back( Range( 1, 100 ) );
	lRanges_2.push_back( Range( 151, 200 ) );
	lRanges_2.push_back( Range( 251, 300 ) );
	RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );

	bool exp = true;
	bool obs = ras_1.isStrictlyContained( ras_2 );
	CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_isStrictlyContained_differentSizes_true_3( void )
{
	RangeAlign ra_1 = RangeAlign();
	std::list<Range> lRanges_1;
	lRanges_1.push_back( Range( 21, 90 ) );
	lRanges_1.push_back( Range( 155, 190 ) );
	lRanges_1.push_back( Range( 160, 195 ) );
	RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );

	RangeAlign ra_2 = RangeAlign();
	std::list<Range> lRanges_2;
	lRanges_2.push_back( Range( 1, 100 ) );
	lRanges_2.push_back( Range( 151, 200 ) );
	RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );

	bool exp = true;
	bool obs = ras_1.isStrictlyContained( ras_2 );
	CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_isStrictlyContained_differentSizes_false( void )
{
	RangeAlign ra_1 = RangeAlign();
	std::list<Range> lRanges_1;
	lRanges_1.push_back( Range( 21, 90 ) );
	lRanges_1.push_back( Range( 155, 190 ) );
	lRanges_1.push_back( Range( 195, 200 ) );
	RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );

	RangeAlign ra_2 = RangeAlign();
	std::list<Range> lRanges_2;
	lRanges_2.push_back( Range( 1, 100 ) );
	lRanges_2.push_back( Range( 151, 200 ) );
	RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );

	bool exp = false;
	bool obs = ras_1.isStrictlyContained( ras_2 );
	CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_isContained( void )
{
	RangeAlign ra_1 = RangeAlign();
	std::list<Range> lRanges_1;
	lRanges_1.push_back( Range( 21, 90 ) );
	lRanges_1.push_back( Range( 155, 190 ) );
	RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );

	RangeAlign ra_2 = RangeAlign();
	std::list<Range> lRanges_2;
	lRanges_2.push_back( Range( 21, 100 ) );
	lRanges_2.push_back( Range( 155, 200 ) );
	RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );

	bool exp = true;
	bool obs = ras_1.isContained( ras_2 );
	CPPUNIT_ASSERT_EQUAL( exp, obs );

	lRanges_2.clear();
	lRanges_2.push_back( Range( 21, 100 ) );
	lRanges_2.push_back( Range( 160, 200 ) );
	ras_2 = RangeAlignSet( &ra_2, lRanges_2 );

	exp = false;
	obs = ras_1.isStrictlyContained( ras_2 );
	CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_doesItContain( void )
{
	RangeAlign ra_1 = RangeAlign();
	std::list<Range> lRanges_1;
	lRanges_1.push_back( Range( 21, 100 ) );
	lRanges_1.push_back( Range( 155, 200 ) );
	RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );

	RangeAlign ra_2 = RangeAlign();
	std::list<Range> lRanges_2;
	lRanges_2.push_back( Range( 21, 90 ) );
	lRanges_2.push_back( Range( 161, 200 ) );
	RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );

	bool exp = true;
	bool obs = ras_1.doesItContain( ras_2 );
	CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_merge_set( void )
{
  RangeAlign ra_1 = RangeAlign();
  std::list<Range> lRanges_1;
  lRanges_1.push_back( Range( 11, 90 ) );
  lRanges_1.push_back( Range( 155, 190 ) );
  lRanges_1.push_back( Range( 21, 100 ) );
  lRanges_1.push_back( Range( 131, 170 ) );
  RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );

  RangeAlign ra_exp = RangeAlign();
  std::list<Range> lRanges_exp;
  lRanges_exp.push_back( Range( 11, 100 ) );
  lRanges_exp.push_back( Range( 131, 190 ) );
  RangeAlignSet ras_exp = RangeAlignSet( &ra_exp, lRanges_exp );

  ras_1.merge_set();

  CPPUNIT_ASSERT_EQUAL( ras_exp, ras_1 );
}

void Test_RangeAlignSet::test_merge_checkCoordinates( void )
{
  RangeAlign ra_1 = RangeAlign();
  std::list<Range> lRanges_1;
  lRanges_1.push_back( Range( 11, 90 ) );
  lRanges_1.push_back( Range( 155, 190 ) );
  RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );

  RangeAlign ra_2 = RangeAlign();
  std::list<Range> lRanges_2;
  lRanges_2.push_back( Range( 21, 100 ) );
  lRanges_2.push_back( Range( 131, 170 ) );
  RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );

  RangeAlign ra_exp = RangeAlign();
  std::list<Range> lRanges_exp;
  lRanges_exp.push_back( Range( 11, 100 ) );
  lRanges_exp.push_back( Range( 131, 190 ) );
  RangeAlignSet ras_exp = RangeAlignSet( &ra_exp, lRanges_exp );

  ras_1.merge( ras_2 );
  CPPUNIT_ASSERT_EQUAL( ras_exp, ras_1 );
}

void Test_RangeAlignSet::test_merge_checkIncluded_false_false( void )
{
  RangeAlign ra_1 = RangeAlign();
  std::list<Range> lRanges_1;
  lRanges_1.push_back( Range( 11, 90 ) );
  lRanges_1.push_back( Range( 155, 190 ) );
  RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );
  ras_1.setIncluded( false );

  RangeAlign ra_2 = RangeAlign();
  std::list<Range> lRanges_2;
  lRanges_2.push_back( Range( 21, 100 ) );
  lRanges_2.push_back( Range( 131, 170 ) );
  RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );
  ras_2.setIncluded( false );

  RangeAlign ra_exp = RangeAlign();
  std::list<Range> lRanges_exp;
  lRanges_exp.push_back( Range( 11, 100 ) );
  lRanges_exp.push_back( Range( 131, 190 ) );
  RangeAlignSet ras_exp = RangeAlignSet( &ra_exp, lRanges_exp );
  ras_exp.setIncluded( false );

  ras_1.merge( ras_2 );
  CPPUNIT_ASSERT_EQUAL( ras_exp, ras_1 );
}

void Test_RangeAlignSet::test_merge_checkIncluded_true_true( void )
{
  RangeAlign ra_1 = RangeAlign();
  std::list<Range> lRanges_1;
  lRanges_1.push_back( Range( 11, 90 ) );
  lRanges_1.push_back( Range( 155, 190 ) );
  RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );
  ras_1.setIncluded( true );

  RangeAlign ra_2 = RangeAlign();
  std::list<Range> lRanges_2;
  lRanges_2.push_back( Range( 21, 100 ) );
  lRanges_2.push_back( Range( 131, 170 ) );
  RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );
  ras_2.setIncluded( true );

  RangeAlign ra_exp = RangeAlign();
  std::list<Range> lRanges_exp;
  lRanges_exp.push_back( Range( 11, 100 ) );
  lRanges_exp.push_back( Range( 131, 190 ) );
  RangeAlignSet ras_exp = RangeAlignSet( &ra_exp, lRanges_exp );
  ras_exp.setIncluded( true );

  ras_1.merge( ras_2 );
  CPPUNIT_ASSERT_EQUAL( ras_exp, ras_1 );
}

void Test_RangeAlignSet::test_merge_checkIncluded_sameCoordinates_sameStrand_true_false( void )
{
  RangeAlign ra_1 = RangeAlign();
  std::list<Range> lRanges_1;
  lRanges_1.push_back( Range( 11, 90 ) );
  lRanges_1.push_back( Range( 155, 190 ) );
  RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );
  ras_1.setIncluded( true );

  RangeAlign ra_2 = RangeAlign();
  std::list<Range> lRanges_2;
  lRanges_2.push_back( Range( 155, 190 ) );
  lRanges_2.push_back( Range( 11, 90 ) );
  RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );
  ras_2.setIncluded( false );

  RangeAlign ra_exp = RangeAlign();
  std::list<Range> lRanges_exp;
  lRanges_exp.push_back( Range( 11, 90 ) );
  lRanges_exp.push_back( Range( 155, 190 ) );
  RangeAlignSet ras_exp = RangeAlignSet( &ra_exp, lRanges_exp );
  ras_exp.setIncluded( true );

  ras_1.merge( ras_2 );
  CPPUNIT_ASSERT_EQUAL( ras_exp, ras_1 );
}

void Test_RangeAlignSet::test_merge_checkIncluded_1false_2true_2notContained( void )
{
  RangeAlign ra_1 = RangeAlign();
  std::list<Range> lRanges_1;
  lRanges_1.push_back( Range( 11, 90 ) );
  lRanges_1.push_back( Range( 155, 190 ) );
  RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );
  ras_1.setIncluded( false );

  RangeAlign ra_2 = RangeAlign();
  std::list<Range> lRanges_2;
  lRanges_2.push_back( Range( 21, 100 ) );
  lRanges_2.push_back( Range( 131, 170 ) );
  RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );
  ras_2.setIncluded( true );

  RangeAlign ra_exp = RangeAlign();
  std::list<Range> lRanges_exp;
  lRanges_exp.push_back( Range( 11, 100 ) );
  lRanges_exp.push_back( Range( 131, 190 ) );
  RangeAlignSet ras_exp = RangeAlignSet( &ra_exp, lRanges_exp );
  ras_exp.setIncluded( false );

  ras_1.merge( ras_2 );
  CPPUNIT_ASSERT_EQUAL( ras_exp.getIncluded(), ras_1.getIncluded() );
  CPPUNIT_ASSERT_EQUAL( ras_exp, ras_1 );
}

void Test_RangeAlignSet::test_merge_checkIncluded_1true_2false_2contained( void )
{
  RangeAlign ra_1 = RangeAlign();
  std::list<Range> lRanges_1;
  lRanges_1.push_back( Range( 11, 100 ) );
  lRanges_1.push_back( Range( 131, 190 ) );
  RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );
  ras_1.setIncluded( true );

  RangeAlign ra_2 = RangeAlign();
  std::list<Range> lRanges_2;
  lRanges_2.push_back( Range( 21, 90 ) );
  lRanges_2.push_back( Range( 155, 170 ) );
  RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );
  ras_2.setIncluded( false );

  RangeAlign ra_exp = RangeAlign();
  std::list<Range> lRanges_exp;
  lRanges_exp.push_back( Range( 11, 100 ) );
  lRanges_exp.push_back( Range( 131, 190 ) );
  RangeAlignSet ras_exp = RangeAlignSet( &ra_exp, lRanges_exp );
  ras_exp.setIncluded( true );

  ras_1.merge( ras_2 );
  CPPUNIT_ASSERT_EQUAL( ras_exp.getIncluded(), ras_1.getIncluded() );
  CPPUNIT_ASSERT_EQUAL( ras_exp, ras_1 );
}

void Test_RangeAlignSet::test_merge_checkIncluded_1true_2false_2notContained( void )
{
  RangeAlign ra_1 = RangeAlign();
  std::list<Range> lRanges_1;
  lRanges_1.push_back( Range( 11, 90 ) );
  lRanges_1.push_back( Range( 155, 190 ) );
  RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );
  ras_1.setIncluded( true );

  RangeAlign ra_2 = RangeAlign();
  std::list<Range> lRanges_2;
  lRanges_2.push_back( Range( 21, 100 ) );
  lRanges_2.push_back( Range( 131, 170 ) );
  RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );
  ras_2.setIncluded( false );

  RangeAlign ra_exp = RangeAlign();
  std::list<Range> lRanges_exp;
  lRanges_exp.push_back( Range( 11, 100 ) );
  lRanges_exp.push_back( Range( 131, 190 ) );
  RangeAlignSet ras_exp = RangeAlignSet( &ra_exp, lRanges_exp );
  ras_exp.setIncluded( false );

  ras_1.merge( ras_2 );
  CPPUNIT_ASSERT_EQUAL( ras_exp.getIncluded(), ras_1.getIncluded() );
  CPPUNIT_ASSERT_EQUAL( ras_exp, ras_1 );
}

void Test_RangeAlignSet::test_merge_checkIncluded_1false_2true_1contained( void )
{
  RangeAlign ra_1 = RangeAlign();
  std::list<Range> lRanges_1;
  lRanges_1.push_back( Range( 21, 90 ) );
  lRanges_1.push_back( Range( 155, 170 ) );
  RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );
  ras_1.setIncluded( false );

  RangeAlign ra_2 = RangeAlign();
  std::list<Range> lRanges_2;
  lRanges_2.push_back( Range( 11, 100 ) );
  lRanges_2.push_back( Range( 131, 190 ) );
  RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );
  ras_2.setIncluded( true );

  RangeAlign ra_exp = RangeAlign();
  std::list<Range> lRanges_exp;
  lRanges_exp.push_back( Range( 11, 100 ) );
  lRanges_exp.push_back( Range( 131, 190 ) );
  RangeAlignSet ras_exp = RangeAlignSet( &ra_exp, lRanges_exp );
  ras_exp.setIncluded( true );

  ras_1.merge( ras_2 );
  CPPUNIT_ASSERT_EQUAL( ras_exp.getIncluded(), ras_1.getIncluded() );
  CPPUNIT_ASSERT_EQUAL( ras_exp, ras_1 );
}

void Test_RangeAlignSet::test_hasSameRanges_true( void )
{
  RangeAlign ra_1 = RangeAlign();
  std::list<Range> lRanges_1;
  lRanges_1.push_back( Range( 1, 100 ) );
  lRanges_1.push_back( Range( 101, 200 ) );
  RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );
  ras_1.sortUsingStrand();
  ras_1.setIncluded( true );

  RangeAlign ra_2 = RangeAlign();
  std::list<Range> lRanges_2;
  lRanges_2.push_back( Range( 1, 100 ) );
  lRanges_2.push_back( Range( 101, 200 ) );
  RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );
  ras_2.sortUsingStrand();
  ras_2.setIncluded( true );

  bool exp = true;
  bool obs = ras_1.hasSameRanges( ras_2 );
  CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_hasSameRanges_true_unsorted( void )
{
  RangeAlign ra_1 = RangeAlign();
  std::list<Range> lRanges_1;
  lRanges_1.push_back( Range( 1, 100 ) );
  lRanges_1.push_back( Range( 101, 200 ) );
  RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );
  ras_1.sortUsingStrand();
  ras_1.setIncluded( true );

  RangeAlign ra_2 = RangeAlign();
  std::list<Range> lRanges_2;
  lRanges_2.push_back( Range( 101, 200 ) );
  lRanges_2.push_back( Range( 1, 100 ) );
  RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );
  ras_2.sortUsingStrand();
  ras_2.setIncluded( true );

  bool exp = true;
  bool obs = ras_1.hasSameRanges( ras_2 );
  CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_hasSameRanges_true_diffIncluded( void )
{
  RangeAlign ra_1 = RangeAlign();
  std::list<Range> lRanges_1;
  lRanges_1.push_back( Range( 1, 100 ) );
  lRanges_1.push_back( Range( 101, 200 ) );
  RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );
  ras_1.sortUsingStrand();
  ras_1.setIncluded( true );

  RangeAlign ra_2 = RangeAlign();
  std::list<Range> lRanges_2;
  lRanges_2.push_back( Range( 1, 100 ) );
  lRanges_2.push_back( Range( 101, 200 ) );
  RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );
  ras_2.sortUsingStrand();
  ras_2.setIncluded( false );

  bool exp = true;
  bool obs = ras_1.hasSameRanges( ras_2 );
  CPPUNIT_ASSERT_EQUAL( exp, obs );
}

void Test_RangeAlignSet::test_hasSameRanges_falseCoords( void )
{
  RangeAlign ra_1 = RangeAlign();
  std::list<Range> lRanges_1;
  lRanges_1.push_back( Range( 1, 100 ) );
  lRanges_1.push_back( Range( 101, 200 ) );
  RangeAlignSet ras_1 = RangeAlignSet( &ra_1, lRanges_1 );
  ras_1.sortUsingStrand();
  ras_1.setIncluded( true );

  RangeAlign ra_2 = RangeAlign();
  std::list<Range> lRanges_2;
  lRanges_2.push_back( Range( 1, 100 ) );
  lRanges_2.push_back( Range( 110, 200 ) );
  RangeAlignSet ras_2 = RangeAlignSet( &ra_2, lRanges_2 );
  ras_2.sortUsingStrand();
  ras_2.setIncluded( true );

  bool exp = false;
  bool obs = ras_1.hasSameRanges( ras_2 );
  CPPUNIT_ASSERT_EQUAL( exp, obs );
}
