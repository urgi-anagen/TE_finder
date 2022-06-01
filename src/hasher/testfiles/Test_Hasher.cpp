#include <SDGMemBioSeq.h>

#include "Test_Hasher.h"

CPPUNIT_TEST_SUITE_REGISTRATION(Test_Hasher);

//-------------------------------------------------------------------------
void
Test_Hasher::run_test_search_wSW(unsigned int alg, unsigned int verbosity, unsigned int start, unsigned int end, unsigned int numseq,
                                 unsigned int connect_dist, unsigned int min_frag_size, unsigned int min_count,
                                 std::list<RangePair> &frag_list, unsigned int kmer_size, unsigned int mask_hole_period,
                                 unsigned int mask_hole_length, unsigned int kmer_window, unsigned int kmer_dist,
                                 unsigned int bkmer_size, unsigned int step_q, double count_cutoff, double diversity_cutoff,
                                 double gap_pen, bool valid_idx_file) {

    Hasher hsrch(kmer_size, mask_hole_period, mask_hole_length, kmer_window, bkmer_size, kmer_dist,
                 0, min_frag_size, step_q, gap_pen, alg);
    hsrch.load("subject_test.fa", kmer_size, mask_hole_period, mask_hole_length,kmer_window, bkmer_size, kmer_size / 2,
               count_cutoff, diversity_cutoff,
               min_count,valid_idx_file, true);

    end= this->seq.size() - 1;
    hsrch.search(this->seq, start, end, numseq, connect_dist,
                 min_frag_size, false, frag_list, verbosity);

    if(verbosity>0) {
        for (const auto &it : frag_list) {
            it.view();
        }
    }

    // run_test_search_wSW on complement
    BioSeq comp_seq= this->seq.complement();
    std::list< RangePair > compfrag_list;
    hsrch.search(comp_seq, start, end, numseq, connect_dist,
                 min_frag_size, false, compfrag_list, verbosity);


    for (auto rp : compfrag_list) {
        rp.getRangeQ().translate_comp(comp_seq.length());
        if (verbosity > 0) rp.view();
        frag_list.push_back(rp);
    }

    double pen_join=0.01;
    hsrch.fragSeqAlign(frag_list,"query_test.fa","subject_test.fa",false,verbosity);

    if(pen_join>0.0){
        std::cout << "Join fragment with penality " << pen_join << " ..." << std::flush;
        hsrch.fragJoin(frag_list);
        std::cout<<" done!"<<std::endl;
    }
    unsigned ext_len=20;
    hsrch.fragSeqSWAlign(frag_list, ext_len, "query_test.fa", "subject_test.fa", false, verbosity);
}
//-------------------------------------------------------------------------
void
Test_Hasher::run_test_search_wExt(unsigned int alg, unsigned int verbosity, unsigned int start, unsigned int end, unsigned int numseq,
                                 unsigned int connect_dist, unsigned int min_frag_size, unsigned int min_count,
                                 std::list<RangePair> &frag_list, unsigned int kmer_size, unsigned int mask_hole_period,
                                 unsigned int mask_hole_length, unsigned int kmer_window, unsigned int kmer_dist,
                                 unsigned int bkmer_size, unsigned int step_q, double count_cutoff, double diversity_cutoff,
                                 double gap_pen, bool valid_idx_file) {

    Hasher hsrch(kmer_size, mask_hole_period, mask_hole_length, kmer_window, bkmer_size, kmer_dist,
                 0, min_frag_size, step_q, gap_pen, alg);
    hsrch.load("subject_test.fa", kmer_size, mask_hole_period, mask_hole_length,kmer_window, bkmer_size, kmer_size / 2,
               count_cutoff, diversity_cutoff,
               min_count,valid_idx_file, true);

    end= this->seq.size() - 1;
    hsrch.search(this->seq, start, end, numseq, connect_dist,
                 min_frag_size, false, frag_list, verbosity);

    if(verbosity>0) {
        for (const auto &it : frag_list) {
            it.view();
        }
    }

    // run_test_search_wSW on complement
    BioSeq comp_seq= this->seq.complement();
    std::list< RangePair > compfrag_list;
    hsrch.search(comp_seq, start, end, numseq, connect_dist,
                 min_frag_size, false, compfrag_list, verbosity);


    for (auto rp : compfrag_list) {
        rp.getRangeQ().translate_comp(comp_seq.length());
        if (verbosity > 0) rp.view();
        frag_list.push_back(rp);
    }

    double pen_join=0.01;
    hsrch.fragSeqAlign(frag_list,"query_test.fa","subject_test.fa",false,verbosity);

    if(pen_join>0.0){
        std::cout << "Join fragment with penality " << pen_join << " ..." << std::flush;
        hsrch.fragJoin(frag_list);
        std::cout<<" done!"<<std::endl;
    }
    unsigned ext_len=20;
    hsrch.fragSeqExtAlign(frag_list, ext_len, "query_test.fa", "subject_test.fa", false, verbosity);
}
//-------------------------------------------------------------------------
void Test_Hasher::test_search_wSW(void ){
    unsigned alg=0;
    unsigned verbosity=1;

    unsigned start=0,end=0,numseq=1,connect_dist=20,min_frag_size=25,min_count=0;
    std::list< RangePair > frag_list;
    unsigned kmer_size=10, mask_hole_period=0, mask_hole_length=0, kmer_window=20, kmer_dist=1, bkmer_size=2, step_q=1;
    double count_cutoff=1.0, diversity_cutoff=0.0, gap_pen=0.01;
    bool valid_idx_file = false;

    run_test_search_wSW(alg, verbosity, start, end, numseq,
                        connect_dist, min_frag_size, min_count,
                        frag_list, kmer_size, mask_hole_period,
                        mask_hole_length, kmer_window, kmer_dist,
                        bkmer_size, step_q, count_cutoff, diversity_cutoff,
                        gap_pen, valid_idx_file) ;

    std::ostringstream ostr_obs;
    for(auto f: frag_list){
        ostr_obs<<f.getIdentity()<<std::endl;
    }

    std::ostringstream ostr_exp;
    std::string exp="100\n100\n60\n100\n43.9024\n100\n100\n63.4146\n100\n";
    ostr_exp<<exp;

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}
//-------------------------------------------------------------------------
void Test_Hasher::test_searchWHole_wSW(void ){
    unsigned alg=1;
    unsigned verbosity=1;

    unsigned start=0,end=0,numseq=1,connect_dist=20,min_frag_size=25,min_count=0;
    std::list< RangePair > frag_list;
    unsigned kmer_size=10, mask_hole_period=0, mask_hole_length=1, kmer_window=0, kmer_dist=1, bkmer_size=2, step_q=1;
    double count_cutoff=1.0, diversity_cutoff=0.0, gap_pen=0.01;
    bool valid_idx_file = false;


    run_test_search_wSW(alg, verbosity, start, end, numseq,
                        connect_dist, min_frag_size, min_count,
                        frag_list, kmer_size, mask_hole_period,
                        mask_hole_length, kmer_window, kmer_dist,
                        bkmer_size, step_q, count_cutoff, diversity_cutoff,
                        gap_pen, valid_idx_file) ;

    std::ostringstream ostr_obs;
    for(auto f: frag_list){
        ostr_obs<<f.getIdentity()<<std::endl;
    }

    std::ostringstream ostr_exp;
    std::string exp="100\n100\n60\n100\n43.9024\n100\n100\n63.4146\n100\n";
    ostr_exp<<exp;

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());

}
//-------------------------------------------------------------------------
void Test_Hasher::test_searchMinimizer_wSW(void ){

    unsigned verbosity=1;

    unsigned start=0,end=0,numseq=1,connect_dist=20,min_frag_size=25,min_count=0;
    std::list< RangePair > frag_list;
    unsigned kmer_size=10, mask_hole_period=0, mask_hole_length=0, kmer_window=15, kmer_dist=10, bkmer_size=2, step_q=1;
    double count_cutoff=1.0, diversity_cutoff=0.0, gap_pen=0.01;
    bool valid_idx_file = false;
    unsigned alg=2;

    run_test_search_wSW(alg, verbosity, start, end, numseq,
                        connect_dist, min_frag_size, min_count,
                        frag_list, kmer_size, mask_hole_period,
                        mask_hole_length, kmer_window, kmer_dist,
                        bkmer_size, step_q, count_cutoff, diversity_cutoff,
                        gap_pen, valid_idx_file) ;

    std::ostringstream ostr_obs;
    for(auto f: frag_list){
        ostr_obs<<f.getIdentity()<<std::endl;
    }

    std::ostringstream ostr_exp;
    std::string exp="100\n100\n100\n100\n100\n100\n";
    ostr_exp<<exp;

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}
//-------------------------------------------------------------------------
void Test_Hasher::test_search_wExt(void ){
    unsigned alg=0;
    unsigned verbosity=1;

    unsigned start=0,end=0,numseq=1,connect_dist=20,min_frag_size=25,min_count=0;
    std::list< RangePair > frag_list;
    unsigned kmer_size=10, mask_hole_period=0, mask_hole_length=0, kmer_window=20, kmer_dist=1, bkmer_size=2, step_q=1;
    double count_cutoff=1.0, diversity_cutoff=0.0, gap_pen=0.01;
    bool valid_idx_file = false;

    run_test_search_wExt(alg, verbosity, start, end, numseq,
                        connect_dist, min_frag_size, min_count,
                        frag_list, kmer_size, mask_hole_period,
                        mask_hole_length, kmer_window, kmer_dist,
                        bkmer_size, step_q, count_cutoff, diversity_cutoff,
                        gap_pen, valid_idx_file) ;

    std::ostringstream ostr_obs;
    for(auto f: frag_list){
        ostr_obs<<f.getIdentity()<<std::endl;
    }

    std::ostringstream ostr_exp;
    std::string exp="100\n100\n100\n100\n100\n100\n100\n100\n100\n";
    ostr_exp<<exp;

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}
//-------------------------------------------------------------------------
void Test_Hasher::test_searchWHole_wExt(void ){
    unsigned alg=1;
    unsigned verbosity=1;

    unsigned start=0,end=0,numseq=1,connect_dist=20,min_frag_size=25,min_count=0;
    std::list< RangePair > frag_list;
    unsigned kmer_size=10, mask_hole_period=0, mask_hole_length=1, kmer_window=0, kmer_dist=1, bkmer_size=2, step_q=1;
    double count_cutoff=1.0, diversity_cutoff=0.0, gap_pen=0.01;
    bool valid_idx_file = false;


    run_test_search_wExt(alg, verbosity, start, end, numseq,
                        connect_dist, min_frag_size, min_count,
                        frag_list, kmer_size, mask_hole_period,
                        mask_hole_length, kmer_window, kmer_dist,
                        bkmer_size, step_q, count_cutoff, diversity_cutoff,
                        gap_pen, valid_idx_file) ;

    std::ostringstream ostr_obs;
    for(auto f: frag_list){
        ostr_obs<<f.getIdentity()<<std::endl;
    }

    std::ostringstream ostr_exp;
    std::string exp="100\n100\n100\n100\n100\n100\n100\n100\n100\n";
    ostr_exp<<exp;

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());

}
//-------------------------------------------------------------------------
void Test_Hasher::test_searchMinimizer_wExt(void ){

    unsigned verbosity=1;

    unsigned start=0,end=0,numseq=1,connect_dist=20,min_frag_size=25,min_count=0;
    std::list< RangePair > frag_list;
    unsigned kmer_size=10, mask_hole_period=0, mask_hole_length=0, kmer_window=15, kmer_dist=10, bkmer_size=2, step_q=1;
    double count_cutoff=1.0, diversity_cutoff=0.0, gap_pen=0.01;
    bool valid_idx_file = false;
    unsigned alg=2;

    run_test_search_wExt(alg, verbosity, start, end, numseq,
                        connect_dist, min_frag_size, min_count,
                        frag_list, kmer_size, mask_hole_period,
                        mask_hole_length, kmer_window, kmer_dist,
                        bkmer_size, step_q, count_cutoff, diversity_cutoff,
                        gap_pen, valid_idx_file) ;

    std::ostringstream ostr_obs;
    for(auto f: frag_list){
        ostr_obs<<f.getIdentity()<<std::endl;
    }

    std::ostringstream ostr_exp;
    std::string exp="100\n100\n100\n100\n100\n100\n";
    ostr_exp<<exp;

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}
//------------------------------------------------------------------------------------------------------------
void Test_Hasher::test_diagSearchDist( void )
{
    unsigned word_len=10;
    unsigned word_dist=1;
    Hasher hsrch(word_len, word_dist, 1);

    std::vector<std::list< HashDNASeq::Diag > > diag_map(3);
    //HashDNASeq::Diag(diag,pos,seq)
    diag_map[1].push_back(HashDNASeq::Diag(1, 10, 1));
    diag_map[1].push_back(HashDNASeq::Diag(1, 20, 1));
    diag_map[1].push_back(HashDNASeq::Diag(1, 30, 1));
    diag_map[1].push_back(HashDNASeq::Diag(1, 30, 1));
    diag_map[1].push_back(HashDNASeq::Diag(1, 60, 1));

    diag_map[2].push_back(HashDNASeq::Diag(1, 70, 2));
    diag_map[2].push_back(HashDNASeq::Diag(1, 100, 2));
    diag_map[2].push_back(HashDNASeq::Diag(1, 130, 2));
    diag_map[2].push_back(HashDNASeq::Diag(1, 140, 2));

    diag_map[1].push_back(HashDNASeq::Diag(2, 100, 1));
    diag_map[1].push_back(HashDNASeq::Diag(2, 110, 1));
    diag_map[1].push_back(HashDNASeq::Diag(2, 120, 1));

    std::list<RangePair> frag;
    unsigned numseq=1, dist=20, min_frag_len=0,verbose=0;
    hsrch.diagSearchDist(numseq, diag_map,dist,word_len,min_frag_len, frag,verbose);

    frag.sort();
    std::ostringstream ostr_obs;
    for(auto it : frag)
    {
        ostr_obs<<it.getRangeQ().getStart()<<".."<<it.getRangeQ().getEnd()
        <<" "<<it.getRangeS().getStart()<<".."<<it.getRangeS().getEnd()<<std::endl;
    }

    std::list<RangePair> frag_exp;
    // (diag+start+1,diag+end+word_size)
    frag_exp.push_back(RangePair(RangeAlign(1,12,41),RangeAlign(0, 11,40)));
    frag_exp.push_back(RangePair(RangeAlign(2, 132,151),RangeAlign(0, 131,150)));
    frag_exp.push_back(RangePair(RangeAlign(1, 103,132),RangeAlign(0, 101,130)));

    frag_exp.sort();
    std::ostringstream ostr_exp;
    for(auto it : frag_exp)
    {
        ostr_exp<<it.getRangeQ().getStart()<<".."<<it.getRangeQ().getEnd()
        <<" "<<it.getRangeS().getStart()<<".."<<it.getRangeS().getEnd()<<std::endl;
    }

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}
//------------------------------------------------------------------------------------------------------------
void Test_Hasher::test_diagSearchScore( void )
{
    unsigned w=10, msk=100, mask_hole_length=1,  bw=2, wd=1,fd=1, minsize=20,step=5,alg=1;
    Hasher hsrch(w, msk, mask_hole_length, bw, wd,fd, minsize,step,alg);


    std::vector<std::list< HashDNASeq::Diag > > diag_map(3);
    //HashDNASeq::Diag(diag,pos,seq)
    diag_map[1].push_back(HashDNASeq::Diag(1, 10, 1));
    diag_map[1].push_back(HashDNASeq::Diag(1, 20, 1));
    diag_map[1].push_back(HashDNASeq::Diag(1, 30, 1));
    diag_map[1].push_back(HashDNASeq::Diag(1, 60, 1));

    diag_map[2].push_back(HashDNASeq::Diag(1, 70, 2));
    diag_map[2].push_back(HashDNASeq::Diag(1, 100, 2));
    diag_map[2].push_back(HashDNASeq::Diag(1, 130, 2));
    diag_map[2].push_back(HashDNASeq::Diag(1, 140, 2));

    diag_map[1].push_back(HashDNASeq::Diag(2, 100, 1));
    diag_map[1].push_back(HashDNASeq::Diag(2, 110, 1));
    diag_map[1].push_back(HashDNASeq::Diag(2, 120, 1));

    std::list<RangePair> frag;
    unsigned numseq=1, min_frag_len=0,verbose=0;
    hsrch.diagSearchScore(numseq, diag_map,min_frag_len, frag,verbose);

    frag.sort();
    std::ostringstream ostr_obs;
    for(auto it : frag)
    {
        ostr_obs<<it.getRangeQ().getNumChr()<<"/"<<it.getRangeQ().getStart()<<".."<<it.getRangeQ().getEnd()
                <<" "<<it.getRangeS().getNumChr()<<"/"<<it.getRangeS().getStart()<<".."<<it.getRangeS().getEnd()
                <<" score="<<it.getScore()<<std::endl;
    }

    std::list<RangePair> frag_exp;
    // (diag+start+1,diag+end+word_size)
    frag_exp.push_back(RangePair( 1,12,71,1,11,70,20,0,0,0));
    frag_exp.push_back(RangePair( 1,132,151,2,131,150,20,0,0,0));
    frag_exp.push_back(RangePair( 1,103,132,1,101,130,30,0,0,0));

    frag_exp.sort();
    std::ostringstream ostr_exp;
    for(auto it : frag_exp)
    {
        ostr_exp<<it.getRangeQ().getNumChr()<<"/"<<it.getRangeQ().getStart()<<".."<<it.getRangeQ().getEnd()
                <<" "<<it.getRangeS().getNumChr()<<"/"<<it.getRangeS().getStart()<<".."<<it.getRangeS().getEnd()
                <<" score="<<it.getScore()<<std::endl;
    }

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}
//-------------------------------------------------------------------------
void Test_Hasher::test_fragJoin(void)
{
    unsigned join_dist=1;
    Hasher hsrch(10, 100, 1, 2, 1, 1,20,1, join_dist, 1);
    std::list< RangePair > frag;

    frag.push_back(RangePair(1, 10, 50,1,110,150,40,0.0, 0.0));
    frag.push_back(RangePair(1, 60, 150,1,160,200,40,0.0, 0.0));
    frag.push_back(RangePair(1, 10, 40,1,120,140,40,0.0, 0.0));
    frag.push_back(RangePair(1, 10, 50,2,110,150,40,0.0, 0.0));
    frag.push_back(RangePair(1, 60, 150,2,160,200,40,0.0, 0.0));
    frag.push_back(RangePair(1, 10, 40,2,120,140,40,0.0, 0.0));

    hsrch.fragJoin(frag);

    frag.sort(RangePair::less);
    std::ostringstream ostr_obs;
    for(std::list< RangePair >::iterator rp_it=frag.begin(); rp_it!=frag.end();rp_it++)
    {
        rp_it->write_raw(ostr_obs);
    }

    //std::cout<<"\n"<<ostr_obs.str()<<std::endl;
    std::list< RangePair > frag_exp;

    frag_exp.push_back(RangePair(1, 10, 150,1,110,200,70,0.0, 0.0));
    frag_exp.push_back(RangePair(1, 10, 40,1,120,140,40,0.0, 0.0));
    frag_exp.push_back(RangePair(1, 10, 150,2,110,200,70,0.0, 0.0));
    frag_exp.push_back(RangePair(1, 10, 40,2,120,140,40,0.0, 0.0));

    frag_exp.sort(RangePair::less);
    std::ostringstream ostr_exp;
    for(std::list< RangePair >::iterator rp_it=frag_exp.begin(); rp_it!=frag_exp.end();rp_it++)
    {
        rp_it->write_raw(ostr_exp);
    }

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}