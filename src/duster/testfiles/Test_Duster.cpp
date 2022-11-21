#include <SDGMemBioSeq.h>

#include "Test_Duster.h"

CPPUNIT_TEST_SUITE_REGISTRATION(Test_Duster);
//-------------------------------------------------------------------------
void Test_Duster::test_search(void ){

    std::ostringstream ostr;
    ostr<<"ATATTTATTTTAGCGTTTACGCTATGTGTTGCGTATTGCTAATCGCTATGATTATATTTATTTTAGCGTTTACGCTATG";
    ostr<<"TTACGCTATGTGTTATTTTTAGCGTTATTGCTAGCGTTTGCGATATTTATTTAATCGCTATGATTATATTTACGCTATG";
    ostr<<"ATATTTCGCGCTATGTGTTGCGATAGCGTTTATTATACCTATATCGCTATGATTATATTTATTTTTAGCGTTTTGTATG";
    BioSeq seq=ostr.str();
    std::ofstream fout_query("query_test.fa");
    fout_query << ">query_test"<<std::endl<<ostr.str();
    fout_query.close();

    BioSeq subseq1=seq.subseq(10-1,51);
    subseq1.header="test1 10..60";
    BioSeq subseq3=seq.subseq(110-1,51);
    subseq3.header="test1 110..160";
    BioSeq subseq2=seq.subseq(50-1,51);
    subseq2.header="test2 comp 100..50";
    subseq2=subseq2.complement();

    std::ostringstream str_fasta;
    str_fasta << ">" << subseq1.header << std::endl;
    str_fasta << subseq1 << std::endl;
    str_fasta << ">" << subseq2.header << std::endl;
    str_fasta << subseq2 << std::endl;
    str_fasta << ">" << subseq3.header << std::endl;
    str_fasta << subseq3 << std::endl;


    std::ofstream fout_subject("subject_test.fa");
    fout_subject << str_fasta.str();
    fout_subject.close();

    unsigned start=5,end=200,numseq=1,min_frag_size=35,verbosity=0,min_count=0;
    std::vector< std::pair<unsigned,unsigned> > frag,frag_comp,fmerged;
    unsigned kmer_size=10, kmask=12, mask_hole_length=1, kmer_window=0, kmer_dist=1, bkmer_size=2, step_q=1;
    double count_cutoff=1.0, diversity_cutoff=0.0;
    bool valid_idx_file = true;

    Duster dstr(kmer_size, kmask, mask_hole_length, bkmer_size, kmer_dist, 0, min_frag_size, step_q);
    dstr.load("subject_test.fa", kmer_size, kmask, mask_hole_length, kmer_window, bkmer_size, kmer_size / 2,
               count_cutoff, diversity_cutoff,
               min_count,valid_idx_file, true, true);

    dstr.search(seq, numseq, start, end, numseq, fmerged, verbosity);

    for (  std::vector< std::pair<unsigned,unsigned> >::iterator it = fmerged.begin();
         it!=fmerged.end(); it++){
        std::cout<<it->first<<".."<<it->second<<std::endl;
    }

    std::system("rm query_test.fa subject_test.fa subject_test.fa.kidx");
}
void Test_Duster::test_fragMerge(void)
{
	unsigned word_len=10;
	unsigned word_dist=1;
	Duster dstr(word_len, word_dist);

	std::vector< std::pair<unsigned,unsigned> > frag;
	frag.push_back(std::pair<unsigned,unsigned>(10,50));
	frag.push_back(std::pair<unsigned,unsigned>(10,50));
	frag.push_back(std::pair<unsigned,unsigned>(100,120));
	frag.push_back(std::pair<unsigned,unsigned>(110,120));
	frag.push_back(std::pair<unsigned,unsigned>(110,130));
	frag.push_back(std::pair<unsigned,unsigned>(110,145));
	frag.push_back(std::pair<unsigned,unsigned>(120,130));
	frag.push_back(std::pair<unsigned,unsigned>(140,150));
	frag.push_back(std::pair<unsigned,unsigned>(200,300));
	frag.push_back(std::pair<unsigned,unsigned>(250,260));
	frag.push_back(std::pair<unsigned,unsigned>(350,360));
	frag.push_back(std::pair<unsigned,unsigned>(360,370));
	frag.push_back(std::pair<unsigned,unsigned>(400,460));
	frag.push_back(std::pair<unsigned,unsigned>(500,600));
	frag.push_back(std::pair<unsigned,unsigned>(550,600));

	std::vector< std::pair<unsigned,unsigned> > fmerged;
	dstr.fragMerge(frag,(unsigned)((word_dist+1)*word_len),fmerged);

    sort(fmerged.begin(),fmerged.end());
	std::ostringstream ostr_obs;
	unsigned size=fmerged.size();
	for(unsigned i=0; i<size; ++i)
		{
			ostr_obs<<fmerged[i].first<<".."<<fmerged[i].second<<std::endl;
	    }

	//std::cout<<"\n"<<ostr_obs.str()<<std::endl;

	std::vector< std::pair<unsigned,unsigned> > fmerged_exp;
	fmerged_exp.push_back(std::pair<unsigned,unsigned>(10,50));
	fmerged_exp.push_back(std::pair<unsigned,unsigned>(100,150));
	fmerged_exp.push_back(std::pair<unsigned,unsigned>(200,300));
	fmerged_exp.push_back(std::pair<unsigned,unsigned>(350,370));
    fmerged_exp.push_back(std::pair<unsigned,unsigned>(400,460));
	fmerged_exp.push_back(std::pair<unsigned,unsigned>(500,600));

    sort(fmerged_exp.begin(),fmerged_exp.end());
    std::ostringstream ostr_exp;
	size=fmerged_exp.size();
	for(unsigned i=0; i<size; ++i)
		{
			ostr_exp<<fmerged_exp[i].first<<".."<<fmerged_exp[i].second<<std::endl;
	    }

	CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}

void Test_Duster::test_runAsScript( void ){
    SDGString inputFileNameGenome = "DmelChr4.fa";
    SDGString inputFileNameTE = "DmelChr4_denovoLibTEs.fa";
    SDGString expFileName = "expDmelChr4.fa.1.duster.bed";

    SDGString prefixFileName = "test_runAsScript";
    SDGString obsFileName = "DmelChr4.fa.1.duster.bed";
    SDGString diff_result = prefixFileName+"result.txt";

    std::ostringstream cmd;
    cmd<<"../duster"<<std::fixed<<std::setprecision(2)<<VERSION;
    cmd<<" -w 15 -k 4 -d 5 -f 100 -S 7 -n 1 "<<inputFileNameGenome<<" "<<inputFileNameTE;
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

        SDGString file=inputFileNameGenome+".1.duster.bed";
        remove(file.c_str());
        file=inputFileNameGenome+".1.duster.bed.fa";
        remove(file.c_str());
        file=inputFileNameTE+".kidx";
        remove(file.c_str());
    }
}