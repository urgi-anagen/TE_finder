#include   <iostream>
#include <cstdlib>
#include <fstream>
#include <unistd.h> //getopt
#include <string>
#include <map>
#include <list>


#include <SDGString.h>

#include "RangeMap.h"

SDGString outfile = "";
bool include = false;
bool overlap = false;
bool del = false;
bool keep = false;
bool diff = false;
bool merge = false;
unsigned extend = 0;

void help(void) {
    std::cerr << "usage: mapOp"
              << " [<options>]" << std::endl
              << " options:" << std::endl
              << "   -h, --help:\n\t this help" << std::endl
              << "   -s, --subject:\n\t the subject map" << std::endl
              << "   -q, --query:\n\t the query map" << std::endl
              << "   -m, --merge:\n\t merge the ranges in a map" << std::endl
              << "   -e, --extend:\n\t extend the ranges by the size given" << std::endl
              << "   -i, --include:\n\t operate on subject when query is included, default: false" << std::endl
              << "   -o, --overlap:\n\t operate on subject when query overlap, default: true" << std::endl
              << "   -k, --keep:\n\t keep subject when query overlap (option -o) or is included (option -i), default: true"
              << std::endl
              << "   -r, --delete:\n\t delete subject when query  overlap (option -o) or is included (option -i), default: false"
              << std::endl
              << "   -d, --diff:\n\t delete subject coordinate according to query overlap, default: false" << std::endl
              << "   -f, --outfile:\n\t output filename, default: <subject map>.[keepInclude,keepOverlap,deleteInclude,deleteOverlap, ...]"
              << std::endl;
};

int main(int argc, char *argv[]) {

    try {

        SDGString s_map_filename;
        SDGString q_map_filename;

        int c;
        while (1) {
// 	static struct option long_options[] =
// 	{
// 	  {"help",no_argument, 0, 'h'},
// 	  {"query",required_argument, 0, 'q'},
// 	  {"subject",required_argument, 0, 's'},
// 	  {"include",no_argument, 0, 'i'},
// 	  {"overlap",no_argument, 0, 'o'},
// 	  {"keep",no_argument, 0, 'k'},
// 	  {"delete",no_argument, 0, 'r'},
// 	  {"diff",no_argument, 0, 'd'},
// 	  {"merge",no_argument, 0, 'm'},
// 	  {"extend",no_argument, 0, 'e'},
// 	  {"outfile",required_argument, 0, 'f'},
// 	  {0, 0, 0, 0}
// 	};
// 	/* `getopt_long' stores the option index here. */
// 	int option_index = 0;

// 	c = getopt_long (argc, argv, "hq:s:iokrdme:f:",
// 			 long_options, &option_index);

            c = getopt(argc, argv, "hq:s:iokrdme:f:");

            /* Detect the end of the options. */
            if (c == -1)
                break;

            switch (c) {
                case 'h': {
                    help();
                    return 0;
                }
                case 'q': {
                    q_map_filename = optarg;
                    break;
                }
                case 's': {
                    s_map_filename = optarg;
                    break;
                }
                case 'm': {
                    merge = true;
                    break;
                }
                case 'e': {
                    extend = atoi(optarg);
                    break;
                }
                case 'i': {
                    include = true;
                    break;
                }
                case 'o': {
                    overlap = true;
                    break;
                }
                case 'r': {
                    del = true;
                    break;
                }
                case 'k': {
                    keep = true;
                    break;
                }
                case 'd': {
                    diff = true;
                    break;
                }
                case 'f': {
                    outfile = optarg;
                    break;
                }
                case '?':
                    help();
                    return 1;
                default:
                    abort();
            }
        }

        /* Print any remaining command line arguments (not options). */
        if (++optind < argc) {
            help();
            std::cout << "non-option ARGV-elements: " << std::endl;
            while (optind < argc)
                std::cout << argv[optind++] << std::endl;
            return 1;
        }

        //load map
        RangeMap query_map, subject_map;
        if (s_map_filename != "") {
            std::cout << "load subject map..." << std::endl << std::flush;
            subject_map.load(s_map_filename);
            subject_map.view();
            std::cout << "\t" << subject_map.getCountRange() << " found" << std::endl << std::flush;
        }
        if (q_map_filename != "") {
            std::cout << "load query map..." << std::endl << std::flush;
            query_map.load(q_map_filename);
            std::cout << "\t" << query_map.getCountRange() << " found" << std::endl << std::flush;
        }

        RangeMap mapout;
        if (keep) {
            if (include) {
                subject_map.selectInclude(query_map, mapout);
                mapout.view();
                std::cout << "save..." << std::endl << std::flush;
                if (outfile == "")
                    mapout.save(SDGString(s_map_filename).afterlast("/")
                                + ".keepInclude");
                else
                    mapout.save(outfile);

            }
            if (overlap) {
                subject_map.selectOverlap(query_map, mapout);
                mapout.view();
                std::cout << "save..." << std::endl << std::flush;
                if (outfile == "")
                    mapout.save(SDGString(s_map_filename).afterlast("/")
                                + ".keepOverlap");
                else
                    mapout.save(outfile);
            }
        }
        if (del) {
            if (include) {
                subject_map.selectInclude(query_map, mapout);
                subject_map.diff(mapout);
                std::cout << "save..." << std::endl << std::flush;
                if (outfile == "")
                    subject_map.save(SDGString(s_map_filename).afterlast("/")
                                     + ".deleteInclude");
                else
                    subject_map.save(outfile);

            }
            if (overlap) {
                subject_map.selectOverlap(query_map, mapout);
                subject_map.diff(mapout);
                std::cout << "save..." << std::endl << std::flush;
                if (outfile == "")
                    subject_map.save(SDGString(s_map_filename).afterlast("/")
                                     + ".deleteOverlap");
                else
                    subject_map.save(outfile);

            }
        }
        if (diff) {
            subject_map.diff(query_map);
            std::cout << "save..." << std::endl << std::flush;
            if (outfile == "")
                subject_map.save(SDGString(s_map_filename).afterlast("/")
                                 + ".diff");
            else
                subject_map.save(outfile);
        }

        if (merge) {
            if (s_map_filename != "") {
                subject_map.merge();
                std::cout << "save..." << std::endl << std::flush;
                if (outfile == "")
                    subject_map.save(SDGString(s_map_filename).afterlast("/")
                                     + ".merge");
                else
                    subject_map.save(outfile);
            }
            if (q_map_filename != "") {
                query_map.merge();
                std::cout << "save..." << std::endl << std::flush;
                if (outfile == "")
                    query_map.save(SDGString(q_map_filename).afterlast("/")
                                   + ".merge");
                else
                    query_map.save(outfile);
            }
        }
        if (extend) {
            if (s_map_filename != "") {
                subject_map.extend(extend);
                std::cout << "save..." << std::endl << std::flush;
                if (outfile == "")
                    subject_map.save(SDGString(s_map_filename).afterlast("/")
                                     + ".extend" + SDGString(extend));
                else
                    subject_map.save(outfile);
            }
            if (q_map_filename != "") {
                query_map.extend(extend);
                std::cout << "save..." << std::endl << std::flush;
                if (outfile == "")
                    query_map.save(SDGString(q_map_filename).afterlast("/")
                                   + ".extend" + SDGString(extend));
                else
                    query_map.save(outfile);
            }
        }
    }

    catch (SDGException e) {
        std::cerr << e.message << std::endl;
    }
    exit(EXIT_SUCCESS);
};
