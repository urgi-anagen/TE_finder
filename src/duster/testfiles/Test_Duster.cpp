#include <SDGMemBioSeq.h>

#include "Test_Duster.h"

CPPUNIT_TEST_SUITE_REGISTRATION(Test_Duster);

void Test_Duster::test_fragMerge(void)
{
	unsigned word_len=10;
	unsigned word_dist=1;
	Duster hsrch(word_len, word_dist);

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
	frag.push_back(std::pair<unsigned,unsigned>(450,460));
	frag.push_back(std::pair<unsigned,unsigned>(500,600));
	frag.push_back(std::pair<unsigned,unsigned>(550,600));

	std::vector< std::pair<unsigned,unsigned> > fmerged;
	hsrch.fragMerge(frag,(unsigned)((word_dist+1)*word_len),fmerged);

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
    cmd<<"../../../cmake-build-debug/src/duster/duster"<<std::fixed<<std::setprecision(2)<<VERSION;
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