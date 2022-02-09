/****
 *
 * \file test_u.cpp
 *
 * \brief Unitary tests for duster
 *
 ***/

#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/BriefTestProgressListener.h>
#include "Test_Hasher.h"

int main(int argc, char *argv[]) {

    try{
        // informs test-listener about testresults
        CPPUNIT_NS::TestResult testresult;

        // register listener for collecting the test-results
        CPPUNIT_NS::TestResultCollector collectedresults;
        testresult.addListener(&collectedresults);

        // register listener for per-test progress output
        CPPUNIT_NS::BriefTestProgressListener progress;
        testresult.addListener(&progress);

        // insert test-suite at test-runner by registry
        CPPUNIT_NS::TestRunner testrunner;
        testrunner.addTest(CPPUNIT_NS::TestFactoryRegistry::getRegistry().makeTest());
        testrunner.run(testresult);

        // output results in compiler-format
        CPPUNIT_NS::CompilerOutputter compileroutputter(&collectedresults, std::cerr);
        compileroutputter.write();

        // return 0 if tests were successful
        return collectedresults.wasSuccessful() ? 0 : 1;
    }
    catch (const SDGException &e) {
        std::cerr << "******Exception catched: " << e.message << " ******" << std::endl;
        exit(EXIT_FAILURE);
    }
    catch (const std::exception &e) {
        std::cout << "Caught exception \"" << e.what() << "\"\n";
    }
    catch (const char *msg) {
        std::cerr << msg << std::endl;
        exit(EXIT_FAILURE);
    }
    catch (...) {
        std::cerr << "****** unknown exception catch !!! ******" << std::endl;
        exit(EXIT_FAILURE);
    }

}

