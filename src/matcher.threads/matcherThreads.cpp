/***
 * 
 * matcher.cpp
 *
 ***/


#include <iostream>
#include <cstdlib>

#include "BLRMatcherThreadsParameter.h"
#include "SDGString.h"
#include "Reference.h"
#include "BLRMatcherThreads.h"

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

        BLRMatcherThreads matcherThreads;
        matcherThreads.process(para);

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
