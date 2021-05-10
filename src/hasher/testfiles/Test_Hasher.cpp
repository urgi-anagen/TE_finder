#include <SDGMemBioSeq.h>

#include "Test_Hasher.h"

CPPUNIT_TEST_SUITE_REGISTRATION(Test_Hasher);

//-------------------------------------------------------------------------
void Test_Hasher::test_search(void ){

    unsigned verbosity=1;

    std::ostringstream ostr;
    ostr<<"ATATTTATTTTAGCGTTTACGCTATGTGTTGCGTATTGCTAATCGCTATGATTATATTTATTTTAGCGTTTACGCTATG";
    ostr<<"TTACGCTATGTGTTATTTTTAGCGTTATTGCTAGCGTTTGCGATATTTATTTAATCGCTATGATTATATTTACGCTATG";
    ostr<<"ATATTTCGCGCTATGTGTTGCGATAGCGTTTATTATACCTATATCGCTATGATTATATTTATTTTTAGCGTTTTGTATG";
    BioSeq seq=BioSeq(ostr.str());
    std::ofstream fout_query("query_test.fa");
    fout_query << ">query_test"<<std::endl<<ostr.str();
    fout_query.close();

    BioSeq subseq1=seq.subseq(10-1,50);
    subseq1.header="test1 10..59";

    BioSeq subseq2=seq.subseq(50-1,50);
    subseq2.header="test2 comp 99..50";
    subseq2=subseq2.complement();

    BioSeq subseq3=seq.subseq(0,50);
    subseq3.header="test3 comp 49..1";
    subseq3=subseq3.complement();

    BioSeq subseq4=seq.subseq(0,50);
    subseq4.header="test4 1..49";

    std::ostringstream str_fasta;
    str_fasta << ">" << subseq1.header << std::endl;
    str_fasta << subseq1 << std::endl;
    str_fasta << ">" << subseq2.header << std::endl;
    str_fasta << subseq2 << std::endl;
    str_fasta << ">" << subseq3.header << std::endl;
    str_fasta << subseq3 << std::endl;
    str_fasta << ">" << subseq4.header << std::endl;
    str_fasta << subseq4 << std::endl;

    std::ofstream fout_subject("subject_test.fa");
    fout_subject << str_fasta.str();
    fout_subject.close();

    unsigned start=0,end,numseq=1,connect_dist=20,min_frag_size=35,min_count=0;
    std::list< RangePair > frag_list;
    unsigned kmer_size=10, mask_hole_period=0, mask_hole_length=1, kmer_dist=1, bkmer_size=2, step_q=1;
    double count_cutoff=1.0, diversity_cutoff=0.0;
    bool valid_idx_file = false;

    Hasher hsrch(kmer_size, mask_hole_period, mask_hole_length, bkmer_size, kmer_dist, 0, min_frag_size, step_q);
    hsrch.load("subject_test.fa", kmer_size, mask_hole_period, mask_hole_length, bkmer_size, kmer_size / 2,
               count_cutoff, diversity_cutoff,
               min_count,valid_idx_file, true);

    end=seq.size();
    hsrch.search(seq, start, end, numseq, connect_dist,
           min_frag_size, false, frag_list, verbosity);

    if(verbosity>0) {
        for (const auto &it : frag_list) {
            it.view();
        }
    }

    // search on complement
    BioSeq comp_seq=seq.complement();
    std::list< RangePair > compfrag_list;
    hsrch.search(comp_seq, start, end, numseq, connect_dist,
                 min_frag_size, false, compfrag_list, verbosity);


    for (auto rp : compfrag_list) {
        rp.getRangeQ().translate_comp(comp_seq.length());
        if (verbosity > 0) rp.view();
        frag_list.push_back(rp);
    }

    if(verbosity>0){
        std::cout<<"\nsubseq1="<<subseq1<<std::endl;
        std::cout<<"subseq2="<<subseq2<<std::endl;
        std::cout<<"subseq3="<<subseq3<<std::endl;
        std::cout<<"subseq4="<<subseq4<<std::endl;
    }
    hsrch.fragSeqAlign(frag_list,"query_test.fa","subject_test.fa",false,verbosity);

    std::ostringstream ostr_obs;
    for(auto f: frag_list){
        ostr_obs<<f.getIdentity()<<std::endl;
    }

    std::ostringstream ostr_exp;
    std::string exp="100\n100\n100\n100\n";
    ostr_exp<<exp;

    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());

    remove("query_test.fa");
    remove("subject_test.fa");
    remove("subject_test.fa.kidx");
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
    frag_exp.push_back(RangePair(RangeAlign(1,12,42),RangeAlign(0, 11,41)));
    frag_exp.push_back(RangePair(RangeAlign(2, 132,152),RangeAlign(0, 131,151)));
    frag_exp.push_back(RangePair(RangeAlign(1, 103,133),RangeAlign(0, 101,131)));

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
    frag_exp.push_back(RangePair( 1,12,72,1,11,71,20,0,0,0));
    frag_exp.push_back(RangePair( 1,132,152,2,131,151,20,0,0,0));
    frag_exp.push_back(RangePair( 1,103,133,1,101,131,30,0,0,0));

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
//void Test_Hasher::test_runAsScript(void ){
//    SDGString inputFileNameGenome = "DmelChr4.fa";
//    SDGString inputFileNameTE = "DmelChr4_denovoLibTEs.fa";
//    SDGString expFileName = "expDmelChr4.fa.hasher.align";
//
//    SDGString prefixFileName = "test_runAsScript";
//    SDGString obsFileName = "DmelChr4.fa.final.hasher.align";
//    SDGString diff_result = prefixFileName+"result.txt";
//
//    std::ostringstream cmd;
//    cmd<<"../hasher"<<std::fixed<<std::setprecision(2)<<VERSION;
//    cmd<<" -w 15 -k 4 -d 5 -S 7 -n 1 -s 50 -c 100 -v 1 "<<inputFileNameGenome<<" "<<inputFileNameTE
//            <<" &> "<<inputFileNameGenome<<".log";
//    std::system(cmd.str().c_str());
//
//    std::ostringstream cmd_diff;
//    cmd_diff<<"diff --side-by-side --suppress-common-lines "<<obsFileName<<" "<<expFileName<<" > "<<diff_result;
//    std::system(cmd_diff.str().c_str());
//
//    std::ostringstream obsStr;
//    std::ifstream fin_obs(diff_result);
//    char buff[2048];
//    while(fin_obs.getline(buff,2047,'\n'))
//        obsStr<<buff<<std::endl;
//
//    bool condition=(obsStr.str()=="");
//    CPPUNIT_ASSERT_MESSAGE("Files "+obsFileName+" and "+expFileName+" are differents",condition);
//    if(condition) {
//        remove(diff_result.c_str());
//        remove(obsFileName.c_str());

//        SDGString file=inputFileNameGenome+".1.duster.bed";
//        remove(file.c_str());
//        file=inputFileNameGenome+".1.duster.bed.fa";
//        remove(file.c_str());
//        file=inputFileNameTE+".kidx";
//        remove(file.c_str());
//    }
//}
////-------------------------------------------------------------------------
//void Test_Hasher::test_fragMerge(void)
//{
//    unsigned word_len=10;
//    unsigned word_dist=1;
//    Hasher hsrch(word_len, word_dist);
//
//    std::list< RangePair > frag;
//
//    frag.push_back(RangePair(1, 10, 50,1,110,150,40,0.0, 0.0));
//    frag.push_back(RangePair(1, 100, 150,1,110,150,40,0.0, 0.0));
//    frag.push_back(RangePair(1, 110, 140,1,120,140,40,0.0, 0.0));
//
//    hsrch.fragMerge(frag);
//
//    frag.sort(RangePair::greater);
//    std::ostringstream ostr_obs;
//    for(std::list< RangePair >::iterator rp_it=frag.begin(); rp_it!=frag.end();rp_it++)
//    {
//        rp_it->writetxt(ostr_obs);
//    }
//
//    //std::cout<<"\n"<<ostr_obs.str()<<std::endl;
//    std::list< RangePair > frag_exp;
//
//    frag_exp.push_back(RangePair(1, 10, 50,1,110,150,40,0.0, 0.0));
//    frag_exp.push_back(RangePair(1, 100, 150,1,110,150,40,0.0, 0.0));
//
//    frag_exp.sort(RangePair::greater);
//    std::ostringstream ostr_exp;
//    for(std::list< RangePair >::iterator rp_it=frag_exp.begin(); rp_it!=frag_exp.end();rp_it++)
//    {
//        rp_it->writetxt(ostr_exp);
//    }
//
//    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
//}