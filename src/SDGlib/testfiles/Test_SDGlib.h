/*
 * Test_SDGlib.h
 *
 *  Created on: 5 juil 2018
 *      Author: hquesnev
 */

#ifndef TEST_SDGLIB_H_
#define TEST_SDGLIB_H_
#include <cppunit/extensions/HelperMacros.h>

#include <list>
#include <vector>

class Test_SDGlib : public CppUnit::TestFixture {

	CPPUNIT_TEST_SUITE(Test_SDGlib);

    CPPUNIT_TEST(test_BioSeq);
    CPPUNIT_TEST(test_BioSeq_subseq);
    CPPUNIT_TEST(test_BioSeq_complement);
    CPPUNIT_TEST(test_BioSeq_reverse);
    CPPUNIT_TEST(test_FastaIstream);
    CPPUNIT_TEST(test_FastaOstream);

	CPPUNIT_TEST( test_SDGMemBioSeq );
	CPPUNIT_TEST( test_SDGMemBioSeq_subseq );
	CPPUNIT_TEST( test_SDGMemBioSeq_complement );
	CPPUNIT_TEST( test_SDGMemBioSeq_reverse );
	CPPUNIT_TEST( test_SDGFastaBioSeq );
	CPPUNIT_TEST( test_SDGFastaIstream );
	CPPUNIT_TEST( test_SDGBioSeqDB );

	CPPUNIT_TEST_SUITE_END();

public:

	Test_SDGlib(void) {}

	void setUp()
	{}
	void tearDown()
	{}

protected:
    void test_BioSeq(void);
    void test_BioSeq_subseq(void);
    void test_BioSeq_complement(void);
    void test_BioSeq_reverse(void);
    void test_FastaIstream( void );
    void test_FastaOstream( void );

	void test_SDGMemBioSeq(void);
    void test_SDGMemBioSeq_subseq(void);
    void test_SDGMemBioSeq_complement(void);
    void test_SDGMemBioSeq_reverse(void);
	void test_SDGFastaBioSeq( void );
	void test_SDGFastaIstream( void );
    void test_SDGBioSeqDB( void );

};


#endif /* TEST_SDGLIB_H_ */
