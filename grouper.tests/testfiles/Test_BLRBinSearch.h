/**
 * Generated CppUnit test driver template.
 * To build it, add the following line at the end of
 * your existing Makefile:
 *    include grouper/Test_BLRBinSearch.mk
 * Build the grouper/Test_BLRBinSearch.mk target from the Make Target view
 */

#ifndef TEST_BLRBinSearch
#define TEST_BLRBinSearch
#include <cppunit/extensions/HelperMacros.h>
#include "../BLRBinSearch.h"
#include "../BLRExhaustSearch.h"
#include <list>


class Test_BLRBinSearch : public CppUnit::TestFixture {
	BLRBinSearch blrBs;
	BLRExhaustSearch blrEs;

	CPPUNIT_TEST_SUITE(Test_BLRBinSearch);
//	CPPUNIT_TEST(testsearchFound);
//	CPPUNIT_TEST(testsearchNotFound);
	CPPUNIT_TEST_SUITE_END();

public:

	Test_BLRBinSearch(void) : blrBs(),blrEs() {}

	void setUp()
	{
		std::list<std::pair<unsigned,unsigned> > lseg;

		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(100,200));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(250,350));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(300,400));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(900,1000));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(1500,2000));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(1800,2200));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(2100,2400));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(15000,20000));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(18000,22000));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(21000,24000));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(150000,200000));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(180000,220000));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(210,240000));

		unsigned count=0;
		for(std::list<std::pair<unsigned,unsigned> >::iterator p=lseg.begin();p!=lseg.end();p++)
		{
			blrEs.insert(p->first,p->second,++count);
			blrBs.insert(p->first,p->second,count);
		}
	}
	void tearDown()
	{
	}
protected:
	void show_vector(const std::vector<unsigned>& v)
	{
		std::cout<<"vector:";
		for(std::vector<unsigned>::const_iterator i=v.begin();i!=v.end();i++)
			std::cout<<*i<<" ";
		std::cout<<std::endl;
	}

	void testsearchFound()
	{
		std::list<std::pair<unsigned,unsigned> > lseg;

		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(150,190));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(150,290));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(150,390));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(950,1100));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(1900,3000));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(5000,10000));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(10,20));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(150000,200000));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(500000000,500220000));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(21000,24000));
		lseg.insert(lseg.begin(),std::pair<unsigned,unsigned>(1500,1500));

		std::vector<unsigned> vb ;
		std::vector<unsigned> ve ;
		for(std::list<std::pair<unsigned,unsigned> >::iterator p=lseg.begin();p!=lseg.end();p++)
		{
			ve=blrEs.search(p->first,p->second);
			vb=blrBs.search(p->first,p->second);
			sort(vb.begin(),vb.end());
			sort(ve.begin(),ve.end());
//			std::cout<<"\nbin search ";
//			show_vector(vb);
//			std::cout<<"exhaustive search ";
//			show_vector(ve);
			CPPUNIT_ASSERT(vb.size()==ve.size());
			CPPUNIT_ASSERT(equal(vb.begin(), vb.end(),ve.begin()));
		}
	}

	void testsearchNotFound()
		{
			std::vector<unsigned> v = blrBs.search(10,12);
			//show_vector(v);
			CPPUNIT_ASSERT(v.empty());
		}


};
CPPUNIT_TEST_SUITE_REGISTRATION(Test_BLRBinSearch);
#endif
