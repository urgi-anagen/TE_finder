/***
 * 
 * matcher.threads.cpp
 *
 ***/


#include <iostream>
#include <cstdlib>
#include <thread>

#include "BLRMatcherThreadsParameter.h"
#include "BLRMatchMap.h"
#include "SDGString.h"
#include "Reference.h"
#include "BLRMatcherThreads.h"

//-------------------------------------------------------------------------------------------------------
void loadFromAlignFile(const BLRMatcherThreadsParameter& para, BLRMatchAlign& match_align)
{
    int verbose=para.getVerbose();


    if(verbose>0) std::cout<<"Load the matches..."<<std::endl<<std::flush;
    match_align.load(para, para.getMatchFileName(), verbose);
    if(verbose>0) std::cout<<"Matches were loaded."<<std::endl;

}
//-------------------------------------------------------------------------------------------------------
std::list< std::list<RangePair> > splitInputData(std::list<RangePair>& rp_list, int nb_set)
{
    std::list< std::list<RangePair> > lrpl;
    rp_list.sort( RangePair::less);

    unsigned count_diff_chr=0;
    unsigned numChr,prevNumChr=0;
    for(std::list<RangePair>::iterator rp_list_it=rp_list.begin();rp_list_it!=rp_list.end();rp_list_it++)
    {
        numChr=rp_list_it->getRangeQ().getNumChr();
        if(numChr!=prevNumChr){
            count_diff_chr++;
        }
        prevNumChr=numChr;
    }
    unsigned size=count_diff_chr;
    unsigned set_size= floor((double)size/nb_set);
    std::cout << "set size=" << set_size << " from " << size << " query sequences" << std::endl;

    prevNumChr=0;
    std::list<RangePair>::iterator rp_list_it=rp_list.begin();
    for(int i=0; i < nb_set; i++)
    {
        std::list<RangePair> rpl;
        std::list<RangePair>::iterator rp_list_start=rp_list_it;
        unsigned count_diff_chr=0;
        std::cout<<"Set"<<i+1<<":";
        while(count_diff_chr <= set_size && rp_list_it != rp_list.end())
        {
            numChr=rp_list_it->getRangeQ().getNumChr();
            if(numChr!=prevNumChr){
                count_diff_chr++;
                std::cout<<rp_list_it->getRangeQ().getNameSeq()<<",";
            }
            rp_list_it++;
            prevNumChr=numChr;
        }// move iterator from several position
        if(rp_list_it!=rp_list.end()){
            std::list<RangePair>::iterator rp_list_it_prev=rp_list_it;
            rp_list_it_prev--;
            rpl.splice(rpl.begin(), rp_list, rp_list_start, rp_list_it_prev);
        }
        else
            rpl.splice(rpl.begin(), rp_list, rp_list_start);
        lrpl.push_back(rpl);
        std::cout<<std::endl;
    }
    lrpl.back().splice(lrpl.back().begin(),rp_list);
    return lrpl;
}
//-------------------------------------------------------------------------------------------------------
void runOnSet(std::list<RangePair>& rp_list,  BLRMatcherThreads& mt, int set_num)
{
    std::cout<<"Run set ..."<<set_num<<std::endl;

    mt.process(rp_list);

    std::cout<<"End set ..."<<set_num<<std::endl;
}
//-------------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    try {
        std::cout << "Beginning Matcher (version " << VERSION << ")" << std::endl << std::flush;
        clock_t start, finish;
        double time_spent;
        start = clock();

        BLRMatcherThreadsParameter para;
        para.parseOptArg(argc, argv);
        if (para.getVerbose() > 0)
            para.view(std::cout);

        BLRMatchMap match_map;

        // Read input data
        BLRMatchAlign match_align;
        if (!para.getLoad_path()) // Load from align file
            loadFromAlignFile(para,match_align);
        else
            std::cerr << "No file :" << para.getLoad_path()  << std::endl << std::flush;

        std::list<RangePair> rp_list = match_align.getRpListFromMatchAlign();
        if(para.getVerbose()>0) std::cout<<"list size="<<rp_list.size()<<std::endl;

        unsigned nb_sets=para.getNbThread();
        if(nb_sets==1 && para.getNbSets()>1) nb_sets=para.getNbSets();
        if(para.getVerbose()>0) std::cout<<"...run on "<<nb_sets<<" sets"<<std::endl<<std::flush;
        if(para.getNbSets()!=1 && para.getNbThread()!=1)
            std::cout<<"......Warning: priority of thread number parameters over set number!"<<std::endl<<std::flush;

        if(para.getVerbose()>0) std::cout<<"Split path in "<<nb_sets<<" sets ..."<<std::endl<<std::flush;
        std::list< std::list<RangePair> > lrpl= splitInputData(rp_list, nb_sets);
        if(para.getVerbose()>0) std::cout<<"End split paths"<<std::endl;


        // Initialize matcherThreads objects

        std::list<BLRMatcherThreads*> lmatchers;
        if(para.getVerbose()>0) std::cout<<"Initialization..."<<std::endl<<std::flush;

        for(unsigned i=0; i<nb_sets; i++)
            lmatchers.push_back(new BLRMatcherThreads(&para, &match_map));

        if(para.getVerbose()>0) std::cout<<"Initialization was done."<<std::endl<<std::flush;


        int set_num=0;
        std::list< std::list<RangePair> >::iterator lrpl_it=lrpl.begin();

        if(para.getNbThread()==1) //Unthreaded
        {
            for(std::list<BLRMatcherThreads*>::iterator lmatchers_it=lmatchers.begin()
                    ; lmatchers_it!=lmatchers.end(); lmatchers_it++)
                runOnSet(std::ref(*lrpl_it++), std::ref(**lmatchers_it), ++set_num);
        }
        else // Threaded
        {
            std::vector<std::thread> threads;
            for(std::list<BLRMatcherThreads*>::iterator lmatchers_it=lmatchers.begin()
                    ; lmatchers_it!=lmatchers.end(); lmatchers_it++)
                threads.push_back(std::thread(runOnSet,std::ref(*lrpl_it++), std::ref(**lmatchers_it), ++set_num));
            for(auto &t : threads) t.join();
        }

        // Get the results
        BLRMatchPath map_path;
        if(para.getVerbose()>0) std::cout<<"Concatenate results..."<<std::endl<<std::flush;
        for(std::list<BLRMatcherThreads*>::iterator lmatchers_it=lmatchers.begin()
                ; lmatchers_it!=lmatchers.end(); lmatchers_it++){
            std::list<RangePairSet> rps_list=(*lmatchers_it)->getRpsListFromMapPath();
            delete (*lmatchers_it);
            map_path.setFromRpsList(para,rps_list,0);
        }
        std::cout<<"Total nb of distinct paths ="<<map_path.getNbDistinctPaths()<<std::endl<<std::flush;

        //Write the results
        std::ostringstream filename;
        if (para.getCleanAfter())
            filename << para.getPrefixFileName() << ".clean_match";
        else
            filename << para.getPrefixFileName() << ".match";

        if (para.getVerbose() > 0)
            std::cout << "Write the results..." << std::endl << std::flush;

        SDGString parafile = filename.str() + ".param";
        para.write(filename.str() + ".param");
        map_path.write(filename.str() + ".path");
        map_path.writeBED(filename.str() + ".bed");
        map_path.writeGFF3(filename.str() + ".gff3");

        if (para.getBank() != "<not set!>" && para.getQuery() != "<not set>")
            map_path.writeMatch(para,filename.str() + ".tab", para.getVerbose());
        if (para.getQuery() != "<not set>")
            map_path.writeSeq(para, filename.str() + ".fa", para.getVerbose());
        if (para.getVerbose() > 0)
            std::cout << "Results were written." << std::endl << std::flush;

        finish = clock();
        time_spent = (double) (finish - start) / CLOCKS_PER_SEC;
        std::cout << "End Matcher in " << time_spent << " second. (version " << VERSION << ")" << std::endl;
        exit(EXIT_SUCCESS);
    }// end try

    catch (SDGException e) {
        std::cout << std::endl << e.message;
    }
    catch (...) {
        std::cout << std::endl << "unknown exception catch";
        exit(EXIT_FAILURE);
    }
}
