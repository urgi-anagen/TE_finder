//
// Created by Hadi Quesneville on 2019-01-22.
//

#include "BLRMatcherThreads.h"

//**************************************************************************
void BLRMatcherThreads::process(BLRMatcherThreadsParameter para) {
    clock_t begin, end;
    double time_spent;

    BLRMatchAlign map_align;
    BLRMatchPath map_path;
    BLRMatchJoin matchJoin(para);

    if (para.getVerbose() > 0)
        std::cout << "**Load the matches..." << std::endl << std::flush;
    map_align.load(para, para.getMatchFileName(), para.getVerbose());
    if (para.getVerbose() > 0)
        std::cout << "...Matches were loaded." << std::endl;

    if (para.getJoin_frag()) {

        if (para.getVerbose() > 0)
            std::cout << "**Connect the fragments..." << std::endl << std::flush;
        begin = clock();
        matchJoin.join(map_align, map_path, para.getVerbose());
        end = clock();
        if (para.getVerbose() > 0) {
            time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
            std::cout << "...Fragments were connected in " << time_spent << " seconds." << std::endl;
        }

        map_align.clear(); // save memory
        //computeScoreWithLength(); ?

        if (para.getMerge()) {
            if (para.getVerbose() > 0)
                std::cout << "**Merge overlaping fragments..." << std::endl << std::flush;
            begin = clock();
            matchJoin.merge(map_path, para.getVerbose());
            end = clock();
            if (para.getVerbose() > 0) {
                time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
                std::cout << "...Fragments were merged in " << time_spent << " seconds." << std::endl;
            }
        }

        if (para.getCleanAfter()) {
            if (para.getVerbose() > 0)
                std::cout << "**Clean conflicting fragments..." << std::endl << std::flush;
            begin = clock();
            matchJoin.clean_conflicts(map_path, para.getVerbose());
            end = clock();
            if (para.getVerbose() > 0) {
                time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
                std::cout << "...Fragments were cleaned in " << time_spent << " seconds." << std::endl;
            }

            if (para.getVerbose() > 0)
                std::cout << "**Split conflicting fragments..." << std::endl << std::flush;
            begin = clock();
            matchJoin.split(map_path, para.getVerbose());
            end = clock();
            if (para.getVerbose() > 0) {
                time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
                std::cout << "...Fragments were splitted in " << time_spent << " seconds." << std::endl;
            }
        }


    } else {
        if (para.getVerbose() > 0)
            std::cout << "**No fragment joins !" << std::endl << std::flush;
        matchJoin.noJoin(map_align, map_path, para.getVerbose());
        map_align.clear(); // save memory
    }

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
    map_path.writeBED(filename.str() + ".bed", "black");

    if (para.getBank() != "<not set!>" && para.getQuery() != "<not set>")
        map_path.writeMatch(filename.str() + ".tab", para.getVerbose());
    if (para.getQuery() != "<not set>")
        map_path.writeSeq(filename.str() + ".fa", para.getVerbose());
    if (para.getVerbose() > 0)
        std::cout << "Results were written." << std::endl << std::flush;
}