/**
 * \file grouper.cpp
 */

#include <iostream>
#include <cstdlib>
#include <thread>
#include "SDGError.h"
#include "Reference.h"
#include "SDGString.h"
#include "BlastMatch.h"
#include "BLRGrouperParameter.h"
#include "BLRGroupThreads.h"
#include "BLRGrouperThreads.h"

//-------------------------------------------------------------------------------------------------------
void loadFromAlignFile(BLRGrouperParameter* para, BLRMatchMap& match_map)
{
	int verbose=para->getVerbose();


	if(verbose>0) std::cout<<"Load the matches..."<<std::endl<<std::flush;
    match_map.loadAlign(para->getMatchFileName(), verbose);
	if(verbose>0) std::cout<<"Matches were loaded."<<std::endl;

	if( para->getJoin_frag() ) //Join fragments
		{
			if(verbose>0) std::cout<<"Connect the fragments..."<<std::endl<<std::flush;
			match_map.mapPath(true, false, false, false, verbose-1);
			//match_map.mapPathWithThreads(true, false, false, false, verbose-1);
			if (verbose>0)
					std::cout<<"Fragments were connected."<<std::endl;
		}
	else //No fragment join
	  match_map.mapPath(false, false, false, false, verbose);
}
//-------------------------------------------------------------------------------------------------------
void loadFromPathFile(BLRGrouperParameter* para, BLRMatchMap& match_map)
{
	int verbose=para->getVerbose();

	if(verbose>0) std::cout<<"Load the matches from path file..."<<std::endl<<std::flush;
	//match_map.loadPath();
	match_map.loadPath(para->getPath_filename(),1);
	if(verbose>0) std::cout<<"Load the matches from path file done!"<<std::endl<<std::flush;
}
//-------------------------------------------------------------------------------------------------------
std::list< std::list<RangePairSet> > splitInputData(std::list<RangePairSet>& rp_list, int nb_split)
{
	std::list< std::list<RangePairSet> > lrpl;
	rp_list.sort( RangePair::greaterLengthIdent );

	int size=rp_list.size();
	int chunk=floor(size/nb_split);
	std::cout<<"chunks="<<chunk<<" from "<<size<<" matches"<<std::endl;
	std::list<RangePairSet>::iterator rp_list_it=rp_list.begin();
	for(int i=0; i<nb_split; i++)
	{
		std::list<RangePairSet> rpl;
		for(int j=0; j<chunk && rp_list_it!=rp_list.end(); j++) rp_list_it++; // move iterator from several position
		if(rp_list_it!=rp_list.end())
			rpl.splice(rpl.begin(), rp_list, rp_list.begin(), rp_list_it);
		else
			rpl.splice(rpl.begin(), rp_list, rp_list.begin());
		lrpl.push_back(rpl);
	}
	lrpl.back().splice(lrpl.back().begin(),rp_list);
	return lrpl;
}
//-------------------------------------------------------------------------------------------------------
void runOnSet(std::list<RangePairSet>& rp_list,  BLRGroup& gr, BLRGrouperParameter& para, int set_num)
{
	std::cout<<"Run set ..."<<set_num<<std::endl;

	BLRGrouper grouper(&para);

	// Build groups via simple-link clustering using percentage of overlap
	// between matches and merge matches if necessary
	grouper.group( rp_list, gr,  para.getVerbose() );


	if(para.getVerbose()>1)
		gr.show_group();

	// Select group containing at least n members
	if( para.getSizefilter() > 1 )
		grouper.group_size_filter( gr, para.getVerbose() );

	// Keep groups containing at least n unique members
	if( para.getIncludeFilter() > 1 )
		grouper.include_filter( gr, para.getVerbose() );

	// Build clusters of groups using length of overlap
	if( para.getGraphfilter() > 0 )
		grouper.cluster(gr, para.getVerbose() );

	std::cout<<"End set ..."<<set_num<<std::endl;
}
//-------------------------------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
	try
	{
		std::cout<<"Beginning Grouper (version "<<VERSION<<")"<<std::endl;
		clock_t begin, end;
                double time_spent;
                begin = clock();

        // Get parameters from command line
		BLRGrouperParameter para;
		para.parseOptArg( argc, argv );
		if(para.getVerbose()>0)
			para.view(std::cout);

		// Read input data
		BLRMatchMap match_map(&para);
		if (!para.getLoad_path()) // Load from align file
			loadFromAlignFile(&para,match_map);
		else
			loadFromPathFile(&para,match_map);

		std::list<RangePairSet> rp_list = match_map.getRpsListFromMapPath();
		if(para.getVerbose()>0) std::cout<<"list size="<<rp_list.size()<<std::endl;

		unsigned nb_sets=para.getNbThread();
		if(nb_sets==1 && para.getNbSets()>1) nb_sets=para.getNbSets();
		if(para.getVerbose()>0) std::cout<<"...run on "<<nb_sets<<" sets"<<std::endl<<std::flush;
		if(para.getNbSets()!=1 && para.getNbThread()!=1)
			std::cout<<"......Warning: priority of thread number parameters over set number!"<<std::endl<<std::flush;

		if(para.getVerbose()>0) std::cout<<"Split path in "<<nb_sets<<" sets ..."<<std::endl<<std::flush;
		std::list< std::list<RangePairSet> > lrpl= splitInputData(rp_list, nb_sets);
		if(para.getVerbose()>0) std::cout<<"End split paths"<<std::endl;

		// Initialize grouper grouping objects

		std::list<BLRGroup*> lgroups;
		if(para.getVerbose()>0) std::cout<<"Initialization..."<<std::endl<<std::flush;


		for(unsigned i=0; i<nb_sets; i++)
			lgroups.push_back(new BLRGroup(&para, &match_map, match_map.getNbQseq()));

		if(para.getVerbose()>0) std::cout<<"Initialization was done."<<std::endl<<std::flush;


		int set_num=0;
		std::list< std::list<RangePairSet> >::iterator lrpl_it=lrpl.begin();

		if(para.getNbThread()==1) //Unthreaded
		{
			for(std::list<BLRGroup*>::iterator lgroups_it=lgroups.begin()
					; lgroups_it!=lgroups.end(); lgroups_it++)
				runOnSet(std::ref(*lrpl_it++), std::ref(**lgroups_it),  std::ref(para),++set_num);
		}
		else // Threaded
		{
			std::vector<std::thread> threads;
			for(std::list<BLRGroup*>::iterator lgroups_it=lgroups.begin()
					; lgroups_it!=lgroups.end(); lgroups_it++)
				threads.push_back(std::thread(runOnSet,std::ref(*lrpl_it++), std::ref(**lgroups_it),  std::ref(para),++set_num));
			for(auto &t : threads) t.join();
		}

		// Merge results
		int c=1;
		std::list<BLRGroup*>::iterator lgroups_it=lgroups.begin();
		lgroups_it++;
		BLRGrouper grouper(&para);
		while(lgroups_it!=lgroups.end())
		{
			std::cout<<"Merge group 1 with group "<<++c<<std::endl;
			grouper.mergeGroupsLists(lgroups.front(),*lgroups_it++,para.getVerbose());
		}

		// Save the results in files
		lgroups.front()->save();

		// clean memory
		for(std::list<BLRGroup*>::iterator lgroups_it=lgroups.begin()
				; lgroups_it!=lgroups.end(); lgroups_it++)
			delete *lgroups_it;

		std::cout<<"End Grouper (version "<<VERSION<<")"<<std::endl;
		end = clock();
        time_spent=(double)(end-begin)/CLOCKS_PER_SEC;
        std::cout<<"Time_spent****: "<<time_spent<<std::endl;

		exit( EXIT_SUCCESS );
	}

	catch( SDGException& e )
	{
		std::cout<<std::endl<<e.message;
		exit( EXIT_FAILURE );
	}
	catch(...)
	{
		std::cout<<std::endl<<"unknown exception catch";
		exit( EXIT_FAILURE );
	}
}
