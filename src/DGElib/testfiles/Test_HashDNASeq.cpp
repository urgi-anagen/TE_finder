#include <SDGMemBioSeq.h>

#include "Test_HashDNASeq.h"

CPPUNIT_TEST_SUITE_REGISTRATION(Test_HashDNASeq);
//------------------------------------------------------------------------------------------------------------
void Test_HashDNASeq::test_hashSeqCount( void )
{
	unsigned word_len=2;
	HashDNASeq hsrch(word_len);

	BioSeq seq=BioSeq("ATATTTATTTTAGCGTTTACGCT");
	std::vector<unsigned> word_count((unsigned)pow(4,word_len),0);
    hsrch.hashSeqCountWHole(seq, word_len, word_count);

	std::ostringstream ostr_obs;
	unsigned size=word_count.size();
	for(unsigned i=0; i<size; i++)
		{
			if(word_count[i]!=0)
				ostr_obs<<"["<<i<<"="<<hsrch.hseq.reverse_hash(i)<<"]="<<word_count[i]<<std::endl;
	    }

	//std::cout<<"\n"<<ostr_obs.str()<<std::endl;


    std::ostringstream ostr_exp;
    ostr_exp
		<<"[1=AC]=1\n"
		<<"[2=AG]=1\n"
		<<"[3=AT]=3\n"
		<<"[6=CG]=2\n"
		<<"[7=CT]=1\n"
		<<"[9=GC]=2\n"
		<<"[11=GT]=1\n"
		<<"[12=TA]=4\n"
		<<"[15=TT]=7"<<std::endl;

 
	CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}
//------------------------------------------------------------------------------------------------------------
void Test_HashDNASeq::test_hashSeqCountwHole( void )
{
    unsigned word_len=3;
    unsigned kmask=2;
    HashDNASeq hsrch(word_len, kmask);

    BioSeq seq=BioSeq("CTCTAT");
    std::vector<unsigned> word_count;
    word_count.resize((unsigned)pow(4,hsrch.getEffectiveKmerSize()));

    hsrch.hashSeqCountWHole(seq, word_len, word_count);

    std::ostringstream ostr_obs;
    unsigned size=word_count.size();
    for(unsigned i=0; i<size; i++)
        {
            if(word_count[i]!=0) 
                ostr_obs<<"["<<hsrch.hseq.reverse_hash(i)<<"]="<<word_count[i]<<std::endl;
        }

    std::ostringstream ostr_exp;
    ostr_exp
        <<"[C-A]=1\n"
        <<"[C-C]=1\n"
        <<"[T-T]=2"
        <<std::endl;


    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}
//------------------------------------------------------------------------------------------------------------
void Test_HashDNASeq::test_hash( void )
{
    unsigned word_len=2;
    HashDNASeq hsrch(word_len);

    unsigned key=hsrch.hseq.hash("TT");

    CPPUNIT_ASSERT_EQUAL(unsigned(15),key);
}
//------------------------------------------------------------------------------------------------------------
void Test_HashDNASeq::test_reverse_hash( void )
{
	unsigned word_len=4;
	HashDNASeq hsrch(word_len);

	unsigned key=hsrch.hseq.hash("TTGC");  
	std::string kmer=hsrch.hseq.reverse_hash(key);
 
	CPPUNIT_ASSERT_EQUAL(std::string("TTGC"),kmer);
}
//------------------------------------------------------------------------------------------------------------
void Test_HashDNASeq::test_reverse_hashwHole( void )
{
    unsigned word_len=10;
    unsigned kmask=2; 
    HashDNASeq hsrch(word_len, kmask);
 
    unsigned key=hsrch.hseq.hash("CGTGAGTGGGG"); 
    unsigned key2=hsrch.hseq.hash("CCTCACTCGCG"); 
    std::string kmer=hsrch.hseq.reverse_hash(key);
    
    CPPUNIT_ASSERT_EQUAL(key,key2);
    CPPUNIT_ASSERT_EQUAL(std::string("C-T-A-T-G-"),kmer);
}
//------------------------------------------------------------------------------------------------------------
void Test_HashDNASeq::test_diagSearchDist( void )
{
	unsigned word_len=10;
	unsigned word_dist=1;
	HashDNASeq hsrch(word_len, word_dist, 1);

    std::vector< std::list<HashDNASeq::Diag> > diag_map(3,std::list<HashDNASeq::Diag>(0));

	//HashDNASeq::Diag(diag,pos,seq)
	diag_map[1].push_back(HashDNASeq::Diag(1, 10, 1));
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
	diag_map[1].push_back(HashDNASeq::Diag(2, 100, 1));
	diag_map[1].push_back(HashDNASeq::Diag(2, 110, 1));
	diag_map[1].push_back(HashDNASeq::Diag(2, 120, 1));
	diag_map[1].push_back(HashDNASeq::Diag(2, 120, 1));



	std::vector< std::pair<unsigned,unsigned> > frag;
    hsrch.diagSearchDist(diag_map, (unsigned) ((word_dist + 1) * word_len), word_len, frag);

    sort(frag.begin(),frag.end());
	std::ostringstream ostr_obs;
	unsigned size=frag.size();
	for(unsigned i=0; i<size; ++i)
		{
			ostr_obs<<frag[i].first<<".."<<frag[i].second<<std::endl;
	    }

	//std::cout<<"\n"<<ostr_obs.str()<<std::endl;

	std::vector< std::pair<unsigned,unsigned> > frag_exp;
	// (diag+start+1,diag+end+word_size)
	frag_exp.push_back(std::pair<unsigned,unsigned>(12,41));
	frag_exp.push_back(std::pair<unsigned,unsigned>(132,151));
	frag_exp.push_back(std::pair<unsigned,unsigned>(103,132));

    sort(frag_exp.begin(),frag_exp.end());
    std::ostringstream ostr_exp;
	size=frag_exp.size();
	for(unsigned i=0; i<size; ++i)
		{
			ostr_exp<<frag_exp[i].first<<".."<<frag_exp[i].second<<std::endl;
	    }

	CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}
//------------------------------------------------------------------------------------------------------------
void Test_HashDNASeq::test_minimizer( void )
{
    unsigned kmer=5;
    unsigned window=15;
    MinimizerFuncDNASeq mini(kmer,window);

    unsigned key1=mini.minimizer("ATAGCTTAGTATATCCGGACGTA");
    unsigned key2=mini.minimizer("TACGTAATCGGACTATATATGTA");

    CPPUNIT_ASSERT_EQUAL(key1,key2);
}