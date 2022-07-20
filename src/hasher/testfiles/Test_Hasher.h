/*
 * Test_Duster.h
 *
 *  Created on: 9 nov. 2015
 *      Author: hquesnev
 */

#ifndef TEST_HASHER_H_
#define TEST_HASHER_H_
#include <cppunit/extensions/HelperMacros.h>

#include "HashDNASeq.h"
#include "Hasher.h"
#include <list>
#include <vector>
#include <iomanip>

class Test_Hasher : public CppUnit::TestFixture {

	CPPUNIT_TEST_SUITE(Test_Hasher);

    CPPUNIT_TEST(test_search_wSW );
	CPPUNIT_TEST(test_searchWHole_wSW );
    CPPUNIT_TEST(test_searchMinimizer_wSW );
    CPPUNIT_TEST(test_searchWHoleMinimizer_wSW );

    CPPUNIT_TEST(test_search_wExt );
    CPPUNIT_TEST(test_searchWHole_wExt );
    CPPUNIT_TEST(test_searchMinimizer_wExt );
    CPPUNIT_TEST(test_searchWHoleMinimizer_wExt );

	CPPUNIT_TEST( test_diagSearchDist );
	//CPPUNIT_TEST( test_diagSearchScore );

    CPPUNIT_TEST( test_fragJoin );

	CPPUNIT_TEST_SUITE_END();

    BioSeq seq;
public:

    Test_Hasher(void) {setUp();}
    ~Test_Hasher(void) {tearDown();}

	void setUp()
	{
        std::ostringstream ostr;
        ostr<<"ATATTTATTTTAGCGTTTACGCTATGTGTTGCGTATTGCTAATCGCTATGATTATATTTATTTTAGCGTTTACGCTATG";
        ostr<<"TTACGCTATGTGTTATTTTTAGCGTTATTGCTAGCGTTTGCGATATTTATTTAATCGCTATGATTATATTTACGCTATG";
        ostr<<"ATATTTCGCGCTATGTGTTGCGATAGCGTTTATTATACCTATATCGCTATGATTATATTTATTTTTAGCGTTTTGTATG";
        seq=BioSeq(ostr.str());
        std::cout<<"\nQuery length:"<<seq.size()<<std::endl;
        std::ofstream fout_query("query_test.fa");
        fout_query << ">query_test"<<std::endl<<ostr.str();
        fout_query.close();

        BioSeq subseq1=seq.subseq(10-1,50);
        subseq1.header="test1 10..60";

        BioSeq subseq2=seq.subseq(50-1,50);
        subseq2.header="test2 comp 99..49";
        subseq2=subseq2.complement();

        BioSeq subseq3=seq.subseq(0,50);
        subseq3.header="test3 comp 50..1";
        subseq3=subseq3.complement();

        BioSeq subseq4=seq.subseq(237-50,50);
        subseq4.header="test4 188..237";

        BioSeq subseq5=seq.subseq(0,50);
        subseq5.header="test5 1..49";

        BioSeq subseq6=seq.subseq(237-50,50);
        subseq6.header="test6 comp 237..188";
        subseq6=subseq6.complement();

        std::ostringstream str_fasta;
        str_fasta << ">" << subseq1.header << std::endl;
        str_fasta << subseq1 << std::endl;
        str_fasta << ">" << subseq2.header << std::endl;
        str_fasta << subseq2 << std::endl;
        str_fasta << ">" << subseq3.header << std::endl;
        str_fasta << subseq3 << std::endl;
        str_fasta << ">" << subseq4.header << std::endl;
        str_fasta << subseq4 << std::endl;
        str_fasta << ">" << subseq5.header << std::endl;
        str_fasta << subseq5 << std::endl;
        str_fasta << ">" << subseq6.header << std::endl;
        str_fasta << subseq6 << std::endl;

        std::ofstream fout_subject("subject_test.fa");
        fout_subject << str_fasta.str();
        fout_subject.close();
	}
	void tearDown()
	{
        remove("query_test.fa");
        remove("subject_test.fa");
        remove("subject_test.fa.kidx");
	}

protected:
    void test_search_wSW(void );
    void test_searchWHole_wSW(void );
    void test_searchMinimizer_wSW(void );
    void test_searchWHoleMinimizer_wSW(void );

    void test_search_wExt(void );
    void test_searchWHole_wExt(void );
    void test_searchMinimizer_wExt(void );
    void test_searchWHoleMinimizer_wExt(void );

    void test_diagSearchDist( void );
    //void test_diagSearchScore( void );
    void test_fragJoin( void );

    void run_test_search_wSW(unsigned int alg, unsigned int verbosity, unsigned int start, unsigned int end, unsigned int numseq,
                             unsigned int connect_dist, unsigned int min_frag_size, unsigned int min_count,
                             std::list<RangePair> &frag_list, unsigned int kmer_size, unsigned int mask_hole_period,
                             unsigned int mask_hole_length, unsigned int kmer_window, unsigned int kmer_dist,
                             unsigned int bkmer_size,
                             unsigned int step_q, double count_cutoff, double diversity_cutoff, double gap_pen, bool valid_idx_file);
    void run_test_search_wExt(unsigned int alg, unsigned int verbosity, unsigned int start, unsigned int end, unsigned int numseq,
                             unsigned int connect_dist, unsigned int min_frag_size, unsigned int min_count,
                             std::list<RangePair> &frag_list, unsigned int kmer_size, unsigned int mask_hole_period,
                             unsigned int mask_hole_length, unsigned int kmer_window, unsigned int kmer_dist,
                             unsigned int bkmer_size,
                             unsigned int step_q, double count_cutoff, double diversity_cutoff, double gap_pen, bool valid_idx_file);
};


#endif /* TEST_DUSTER_H_ */
