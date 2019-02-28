#include <SDGMemBioSeq.h>

#include "Test_Duster.h"

CPPUNIT_TEST_SUITE_REGISTRATION(Test_Duster);
//------------------------------------------------------------------------------------------------------------
void Test_Duster::test_hashSeqCount( void )
{
	unsigned word_len=2;
	Duster hsrch(word_len);

	SDGBioSeq seq=newSDGMemBioSeq("ATATTTATTTTAGCGTTTACGCT");
	std::vector<unsigned> word_count;
	word_count.resize((unsigned)pow(4,word_len)); 

	hsrch.hashSeqCount(seq,word_len,word_count);

	std::ostringstream ostr_obs;
	unsigned size=word_count.size();
	for(unsigned i=0; i<size; ++i)
		{
			if(word_count[i]!=0)
				ostr_obs<<"["<<i<<"]="<<word_count[i]<<std::endl;
	    }

	//std::cout<<"\n"<<ostr_obs.str()<<std::endl;


    std::ostringstream ostr_exp;
    ostr_exp
		<<"[1]=1\n"
		<<"[2]=1\n" 
		<<"[3]=3\n"
		<<"[6]=2\n"
		<<"[7]=1\n"
		<<"[9]=2\n"
		<<"[11]=1\n"
		<<"[12]=4\n"
		<<"[15]=7"<<std::endl;

 
	CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
}
//------------------------------------------------------------------------------------------------------------
void Test_Duster::test_hashSeqCountwHole( void )
{
    unsigned word_len=3;
    unsigned kmask=2;
    Duster hsrch(word_len,kmask);

    SDGBioSeq seq=newSDGMemBioSeq("CTCTAT");
    std::vector<unsigned> word_count;
    word_count.resize((unsigned)pow(4,hsrch.getEffectiveKmerSize()));
 
    hsrch.hashSeqCount(seq,word_len,word_count);

    std::ostringstream ostr_obs;
    unsigned size=word_count.size();
    for(unsigned i=0; i<size; ++i)
        {
            if(word_count[i]!=0) 
                ostr_obs<<"["<<hsrch.hseq.reverse_hash(i)<<"]="<<word_count[i]<<std::endl;
        }

    //std::cout<<"\n"<<ostr_obs.str()<<std::endl;

    /* To display all kmers
    word_len=3;
    kmask=100;
    Duster hsrch2(word_len,kmask);

    std::vector<unsigned> word_count2;
    word_count2.resize((unsigned)pow(4,hsrch2.getEffectiveKmerSize()));

    hsrch2.hashSeqCount(seq,word_len,word_count2);

    std::ostringstream ostr_obs2;
    unsigned size2=word_count2.size();
    for(unsigned i=0; i<size2; ++i)
        {
            if(word_count2[i]!=0)
                ostr_obs2<<"["<<hsrch2.hseq.reverse_hash(i)<<"]="<<word_count2[i]<<std::endl;
        }

    std::cout<<"\n"<<ostr_obs2.str()<<std::endl; 
    */

    std::ostringstream ostr_exp;
    ostr_exp
        <<"[C-A]=1\n"
        <<"[C-C]=1\n"
        <<"[T-T]=2"
        <<std::endl;


    CPPUNIT_ASSERT_EQUAL(ostr_exp.str(),ostr_obs.str());
} 
//------------------------------------------------------------------------------------------------------------
void Test_Duster::test_reverse_hash( void )
{
	unsigned word_len=4;
	Duster hsrch(word_len); 

	unsigned key=hsrch.hseq.hash("TTGC");  
	std::string kmer=hsrch.hseq.reverse_hash(key);
 
	CPPUNIT_ASSERT_EQUAL(std::string("TTGC"),kmer);
}
//------------------------------------------------------------------------------------------------------------
void Test_Duster::test_reverse_hashwHole( void )
{
    unsigned word_len=10;
    unsigned kmask=2; 
    Duster hsrch(word_len,kmask); 
 
    unsigned key=hsrch.hseq.hash("CGTGAGTGGGG"); 
    unsigned key2=hsrch.hseq.hash("CCTCACTCGCG"); 
    std::string kmer=hsrch.hseq.reverse_hash(key);
    
    CPPUNIT_ASSERT_EQUAL(key,key2);
    CPPUNIT_ASSERT_EQUAL(std::string("C-T-A-T-G-"),kmer);
}
//------------------------------------------------------------------------------------------------------------
void Test_Duster::test_diagSearch( void )
{
	unsigned word_len=10;
	unsigned word_dist=1;
	Duster hsrch(word_len,word_dist,1);

	std::vector< Duster::Diag > diag_map;

	//Duster::Diag(diag,pos,seq)
	diag_map.push_back(Duster::Diag(1,10,1));
	diag_map.push_back(Duster::Diag(1,10,1));
	diag_map.push_back(Duster::Diag(1,20,1));
	diag_map.push_back(Duster::Diag(1,30,1));
	diag_map.push_back(Duster::Diag(1,30,1));
	diag_map.push_back(Duster::Diag(1,60,1));

	diag_map.push_back(Duster::Diag(1,70,2));
	diag_map.push_back(Duster::Diag(1,100,2));
	diag_map.push_back(Duster::Diag(1,130,2));
	diag_map.push_back(Duster::Diag(1,140,2));

	diag_map.push_back(Duster::Diag(2,100,1));
	diag_map.push_back(Duster::Diag(2,100,1));
	diag_map.push_back(Duster::Diag(2,110,1));
	diag_map.push_back(Duster::Diag(2,120,1));
	diag_map.push_back(Duster::Diag(2,120,1));



	std::vector< std::pair<unsigned,unsigned> > frag;
	hsrch.diagSearch(diag_map,(unsigned)((word_dist+1)*word_len),word_len,frag);

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
void Test_Duster::test_fragMerge( void )
{
	unsigned word_len=10;
	unsigned word_dist=1;
	Duster hsrch(word_len,word_dist);

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
