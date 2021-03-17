#include <SDGMemBioSeq.h>

#include "Test_Hasher.h"

CPPUNIT_TEST_SUITE_REGISTRATION(Test_Hasher);

//-------------------------------------------------------------------------
void Test_Hasher::test_search(void ){

    std::ostringstream ostr;
    ostr<<"ATATTTATTTTAGCGTTTACGCTATGTGTTGCGTATTGCTAATCGCTATGATTATATTTATTTTAGCGTTTACGCTATG";
    ostr<<"TTACGCTATGTGTTATTTTTAGCGTTATTGCTAGCGTTTGCGATATTTATTTAATCGCTATGATTATATTTACGCTATG";
    ostr<<"ATATTTCGCGCTATGTGTTGCGATAGCGTTTATTATACCTATATCGCTATGATTATATTTATTTTTAGCGTTTTGTATG";
    SDGBioSeq seq=newSDGMemBioSeq(ostr.str());
    std::ofstream fout_query("query_test.fa");
    fout_query << ">query_test"<<std::endl<<ostr.str();
    fout_query.close();

    SDGBioSeq subseq1=seq.subseq(10-1,51);
    subseq1.setDE("test1 10..60");
    std::cout<<"\nsubseq1="<<subseq1.toString()<<std::endl;
    SDGBioSeq subseq2=seq.subseq(50-1,51);
    subseq2.setDE("test2 comp 100..50");
    subseq2=subseq2.complement();

    std::ostringstream str_fasta;
    str_fasta << ">" << subseq1.getDE() << std::endl;
    str_fasta << subseq1.toString() << std::endl;;
    str_fasta << ">" << subseq2.getDE() << std::endl;
    str_fasta << subseq2.toString() << std::endl;;

    std::ofstream fout_subject("subject_test.fa");
    fout_subject << str_fasta.str();
    fout_subject.close();

    unsigned start=5,end=200,numseq=1,connect_dist=20,min_frag_size=35,verbosity=0,min_count=0;
    std::list< RangePair > frag_list;
    unsigned kmer_size=10, kmask=11, mask_hole_length=1, kmer_dist=1, bkmer_size=2, step_q=1;
    double count_cutoff=1.0, diversity_cutoff=0.0;
    bool valid_idx_file = true;

    Hasher hsrch(kmer_size, kmask, mask_hole_length, bkmer_size, kmer_dist, 0, min_frag_size, step_q);
    hsrch.load("subject_test.fa", kmer_size, kmer_size, mask_hole_length, bkmer_size, kmer_size / 2,
               count_cutoff, diversity_cutoff,
               min_count,valid_idx_file, true);

    hsrch.search(seq, start, end, numseq, connect_dist,
           min_frag_size, false, frag_list, verbosity);

    for (const auto & it : frag_list){
      it.view();
    }

    // search on complement
    SDGBioSeq comp_seq=seq.complement();
    std::list< RangePair > compfrag_list;
    hsrch.search(comp_seq, start, end, numseq, connect_dist,
                 min_frag_size, false, compfrag_list, verbosity);

    for (auto rp : compfrag_list){
        rp.getRangeQ().translate_comp(comp_seq.length());
        rp.view();
        frag_list.push_back(rp);
    }
    hsrch.fragSeqAlign(frag_list,"query_test.fa","subject_test.fa",false,verbosity);

    std::system("rm query_test.fa subject_test.fa subject_test.fa.kidx");
}
//-------------------------------------------------------------------------
void Test_Hasher::test_runAsScript(void ){
    SDGString inputFileNameGenome = "DmelChr4.fa";
    SDGString inputFileNameTE = "DmelChr4_denovoLibTEs.fa";
    SDGString expFileName = "expDmelChr4.fa.hasher.align";

    SDGString prefixFileName = "test_runAsScript";
    SDGString obsFileName = "DmelChr4.fa.final.hasher.align";
    SDGString diff_result = prefixFileName+"result.txt";

    std::ostringstream cmd;
    cmd<<"../hasher"<<std::fixed<<std::setprecision(2)<<VERSION;
    cmd<<" -w 15 -k 4 -d 5 -S 7 -n 1 -s 50 -c 100 -v 1 "<<inputFileNameGenome<<" "<<inputFileNameTE
            <<" &> "<<inputFileNameGenome<<".log";
    std::system(cmd.str().c_str());

    std::ostringstream cmd_diff;
    cmd_diff<<"diff --side-by-side --suppress-common-lines "<<obsFileName<<" "<<expFileName<<" > "<<diff_result;
    std::system(cmd_diff.str().c_str());

    std::ostringstream obsStr;
    std::ifstream fin_obs(diff_result);
    char buff[2048];
    while(fin_obs.getline(buff,2047,'\n'))
        obsStr<<buff<<std::endl;

    bool condition=(obsStr.str()=="");
    CPPUNIT_ASSERT_MESSAGE("Files "+obsFileName+" and "+expFileName+" are differents",condition);
    if(condition) {
        remove(diff_result.c_str());
        remove(obsFileName.c_str());

//        SDGString file=inputFileNameGenome+".1.duster.bed";
//        remove(file.c_str());
//        file=inputFileNameGenome+".1.duster.bed.fa";
//        remove(file.c_str());
//        file=inputFileNameTE+".kidx";
//        remove(file.c_str());
    }
}
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