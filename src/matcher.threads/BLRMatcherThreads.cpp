//
// Created by Hadi Quesneville on 2019-01-22.
//

#include "BLRMatcherThreads.h"

//**************************************************************************
void BLRMatcherThreads::process(const std::list<RangePair>& rp_list) {
    clock_t begin, end;
    double time_spent;

    BLRMatchAlign map_align;
    BLRMatchJoin matchJoin(*matcher_parameter);

    if (matcher_parameter->getVerbose() > 0)
        std::cout << "**Load the matches..." << std::endl << std::flush;
    map_align.setFromRpsList(*matcher_parameter, rp_list, matcher_parameter->getVerbose());
    if (matcher_parameter->getVerbose() > 0)
        std::cout << "...Matches were loaded." << std::endl;

    if (matcher_parameter->getJoin_frag()) {

        if (matcher_parameter->getVerbose() > 0)
            std::cout << "**Connect the fragments..." << std::endl << std::flush;
        begin = clock();
        matchJoin.join(map_align, map_path, matcher_parameter->getVerbose());
        end = clock();
        if (matcher_parameter->getVerbose() > 0) {
            time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
            std::cout << "...Fragments were connected in " << time_spent << " seconds." << std::endl;
        }

        map_align.clear(); // save memory
        //....computeScoreWithLengthAndId(); ?

        if (matcher_parameter->getMerge()) {
            if (matcher_parameter->getVerbose() > 0)
                std::cout << "**Merge overlaping fragments..." << std::endl << std::flush;
            begin = clock();
            matchJoin.merge(map_path, matcher_parameter->getVerbose());
            end = clock();
            if (matcher_parameter->getVerbose() > 0) {
                time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
                std::cout << "...Fragments were merged in " << time_spent << " seconds." << std::endl;
            }
        }

        if (matcher_parameter->getCleanAfter()) {
            if (matcher_parameter->getVerbose() > 0)
                std::cout << "**Clean conflicting fragments..." << std::endl << std::flush;
            begin = clock();
            matchJoin.clean_conflicts(map_path, matcher_parameter->getVerbose());
            end = clock();
            if (matcher_parameter->getVerbose() > 0) {
                time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
                std::cout << "...Fragments were cleaned in " << time_spent << " seconds." << std::endl;
            }

            if (matcher_parameter->getVerbose() > 0)
                std::cout << "**Split conflicting fragments..." << std::endl << std::flush;
            begin = clock();
            matchJoin.split(map_path, matcher_parameter->getVerbose());
            end = clock();
            if (matcher_parameter->getVerbose() > 0) {
                time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
                std::cout << "...Fragments were splitted in " << time_spent << " seconds." << std::endl;
            }
        }



    } else {
        if (matcher_parameter->getVerbose() > 0)
            std::cout << "**No fragment joins !" << std::endl << std::flush;
        matchJoin.noJoin(map_align, map_path, matcher_parameter->getVerbose());
        map_align.clear(); // save memory
    }


}