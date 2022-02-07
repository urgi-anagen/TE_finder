#include <getopt.h>
#include <fstream>
#include <cstdlib>
#include <list>

#include "SDGString.h"
#include "SDGFastaIstream.h"
#include "FastaOstream.h"
#include "SDGMemBioSeq.h"

#include "Hasher.h"


double filter_cutoff=0.0, pen_join=0.5;
unsigned kmer_size=15, step_q=15, bkmer_size=1, kmer_dist=5,mask_hole_length=1,connect_dist=1,
min_size=20, min_frag_size, chunk_size_kb=0, min_count=0, kmask=4, verbosity=0, overlap=0, nb_iter=1,algorithm=1;

double count_cutoff=1.0, diversity_cutoff=0.0;
bool repeat=false, stat_only=false;
double p_value=0.05;

SDGString outfilename="";

    
void help(void)
{
    std::cerr << "usage: hasher"
            << " [<options>] <fasta query sequence> [<fasta sequence model>]." << std::endl
            << " options:" << std::endl
            << "   -h, --help:\n\t this help" << std::endl
            << "   -w, --kmer:\n\t kmer length, default: " << kmer_size << std::endl
            << "   -S, --step_q:\n\t step on query sequence, default: " << step_q << std::endl
            << "   -k, --kmask:\n\t period of k-mer hole, default: " << kmask << std::endl
            << "   -l, --len_hole_mask:\n\t kmer mask hole length, default: " << mask_hole_length << std::endl
            << "   -d, --kmer_dist:\n\t max number of kmer between two matching kmer to connect, default: "
            << kmer_dist << std::endl
            << "   -j, --pen_join:\n\t penality to join two matching kmer, default: "
            << pen_join << " (enter a value < 0.0 to skip this)" << std::endl
              << "   -s, --min_size:\n\t min size range to report, default: "
              << min_size << std::endl
              << "   -C, --filter_cutoff:\n\t filter kmer with counts over a percentile (Value [0-1]), default: "
              << count_cutoff << std::endl
              << "   -D, --diversity_cutoff:\n\t filter kmer with diversity measure of kmer size used for background probability (Value [0-1]), default: "
              << diversity_cutoff << std::endl
              << "   -m, --min_count:\n\t filter kmer with counts less than this value, default: "
              << min_count << std::endl
              << "   -b, --background_kmer_size:\n\t kmer size to compute kmer background probability, default: "
              << bkmer_size << std::endl
              << "   -p, --p-value:\n\t p-value, default: " << p_value << std::endl
              << "   -o, --file_out:\n\t filename for output," << std::endl
              << "   -c, --chunk_size:\n\t sequence chunk size in kb, default: None" << std::endl
              << "   -n, --nb_iter:\t number of iteration: " << nb_iter << std::endl
              << "   -a, --analysis:\n\t compute kmer statistics only" << std::endl
              << "   -A, --algorithm:\n\t algorithm number (1: kmers with hole, 2: kmer minimizer ), default:" << algorithm << std::endl
              << "   -v, --verbosity:\n\t verbosity level, default:" << verbosity << std::endl;
};
void show_parameter(SDGString filename1,SDGString filename2)
{
    std::cout << "\nRun with parameters:\n"
              << "Query sequences: " << filename1 << std::endl
              << "Model sequences: " << filename2 << std::endl
              << "   -w, --kmer:\t kmer length: " << kmer_size << std::endl
              << "   -S, --step_q:\t step on query sequence: " << step_q << std::endl
              << "   -k, --kmask:\t period of k-mer hole: " << kmask << std::endl
              << "   -l, --len_hole_mask:\t kmer mask hole length: " << mask_hole_length << std::endl
              << "   -d, --kmer_dist:\t max number of kmer between two matching kmer to connect: " << kmer_dist
              << std::endl
              << "   -j, --pen_join:\t penality to join two matching kmer, default: "
              << pen_join << " (enter a value < 0.0 to skip this)" << std::endl
              << "   -s, --min_size:\t min size range to report: " << min_size << std::endl
              << "   -C, --filter_cutoff:\t filter kmer with counts in the last percentile: " << count_cutoff
              << std::endl
              << "   -D, --diversity_cutoff:\t filter kmer with diversity measure of kmer size used for background probability: "
              << diversity_cutoff << std::endl
              << "   -m, --min_count:\t filter kmer with counts less than this value: " << min_count << std::endl
              << "   -b, --background_kmer_size:\t kmer size to compute kmer background probability: " << bkmer_size
              << std::endl
              << "   -p, --p-value:\t p-value, default: " << p_value << std::endl
              << "   -o, --file_out:\t filename for output:" << outfilename << std::endl
              << "   -c, --chunk_size:\t sequence chunk size in kb: " << chunk_size_kb << std::endl
              << "   -n, --nb_iter:\t number of iteration: " << nb_iter << std::endl
              << "   -A, --algorithm:\t algorithm number (1: kmers with hole, 2: kmer minimizer ) : " << algorithm << std::endl
              << "   -v, --verbosity:\t verbosity level: " << verbosity << std::endl;
};
// search on sequence chunk, reverse sequence, and reverse complement
void search_frag(Hasher& hsrch,
                 const BioSeq& seq, const BioSeq& seq_rev,
                 const BioSeq& seq_comp, const BioSeq& seq_revcomp,
                 unsigned start, unsigned end, unsigned numseq, unsigned connect_dist,
                 unsigned min_frag_size, bool repeat,
                 std::list< RangePair >& frag_list,
                 std::list< RangePair >& rev_frag_list,
                 std::list< RangePair >& compfrag_list,
                 std::list< RangePair >& rev_compfrag_list,
                 unsigned verbose)
{
    hsrch.search(seq, start, end, numseq, connect_dist,
                 min_frag_size, repeat, frag_list, verbosity);
    hsrch.search(seq_rev, start, end, numseq, connect_dist,
                 min_frag_size, repeat, rev_frag_list, verbosity);
    hsrch.search(seq_comp, start, end, numseq, connect_dist,
                 min_frag_size, repeat, compfrag_list, verbosity);
    hsrch.search(seq_revcomp, start, end, numseq, connect_dist,
                 min_frag_size, repeat, rev_compfrag_list, verbosity);
}
//translate reverse complement coordinate
void translate_comp(std::list< RangePair >& frag_list,const std::list< RangePair >& compfrag_list,
                    unsigned len_seq){
    for (std::list<RangePair>::const_iterator it = compfrag_list.cbegin();
         it!=compfrag_list.cend(); it++){
        RangePair rp(*it);
        rp.getRangeQ().translate_comp(len_seq);
        frag_list.push_back(rp);
    }
}

int main(int argc, char *argv[]) {
    try {
        std::cout << "Beginning Hasher (version " << VERSION << ")" << std::endl;
        clock_t begin, end;
        double time_spent;
        begin = clock();

        if (argc == 1) {
            help();
            exit(EXIT_SUCCESS);
        }
        int c;
        while (1) {
            static struct option long_options[] =
                    {
                            {"help",                 no_argument,       0, 'h'},
                            {"kmer",                 required_argument, 0, 'w'},
                            {"step_q",               required_argument, 0, 'S'},
                            {"kmask",                required_argument, 0, 'k'},
                            {"len_hole_mask",        required_argument, 0, 'l'},
                            {"kmer_dist",            required_argument, 0, 'd'},
                            {"pen_join",            required_argument, 0, 'j'},
                            {"min_size",             required_argument, 0, 's'},
                            {"filter_cutoff",        required_argument, 0, 'C'},
                            {"min_count",            required_argument, 0, 'm'},
                            {"diversity_cutoff",     required_argument, 0, 'D'},
                            {"background_kmer_size", required_argument, 0, 'b'},
                            {"p_value",              required_argument, 0, 'p'},
                            {"file_out",             required_argument, 0, 'o'},
                            {"chunk_size",           required_argument, 0, 'c'},
                            {"nb_iter",              required_argument, 0, 'n'},
                            {"analysis",             no_argument,       0, 'a'},
                            {"algorithm",            required_argument, 0, 'A'},
                            {"verbosity",            no_argument,       0, 'v'},
                            {0, 0,                                      0, 0}
                    };
            /* `getopt_long' stores the option index here. */
            int option_index = 0;

            c = getopt_long(argc, argv, "hd:f:w:S:k:l:d:j:s:C:D:m:b:p:o:c:n:aA:v:",
                            long_options, &option_index);

            /* Detect the end of the options. */
            if (c == -1)
                break;

            switch (c) {
                case 'h': {
                    help();
                    return 0;
                }
                case 'w': {
                    kmer_size = atoi(optarg);
                    break;
                }
                case 'S': {
                    step_q = atoi(optarg);
                    break;
                }
                case 'k': {
                    kmask = atoi(optarg);
                    break;
                }
                case 'l': {
                    mask_hole_length = atoi(optarg);
                    break;
                }
                case 'd': {
                    kmer_dist = atoi(optarg);
                    break;
                }
                case 'j': {
                    pen_join = atof(optarg);
                    break;
                }
                case 's': {
                    min_size = atoi(optarg);
                    break;
                }
                case 'C': {
                    count_cutoff = atof(optarg);
                    break;
                }
                case 'D': {
                    diversity_cutoff = atof(optarg);
                    break;
                }
                case 'm': {
                    min_count = atoi(optarg);
                    break;
                }
                case 'b': {
                    bkmer_size = atoi(optarg);
                    break;
                }
                case 'p': {
                    p_value = atof(optarg);
                    break;
                }
                case 'o': {
                    outfilename = optarg;
                    break;
                }
                case 'c': {
                    chunk_size_kb = atoi(optarg);
                    break;
                }
                case 'n': {
                    nb_iter = atoi(optarg);
                    break;
                }
                case 'a': {
                    stat_only = true;
                    break;
                }
                case 'A': {
                    algorithm = atoi(optarg);
                    break;
                }
                case 'v': {
                    verbosity = atoi(optarg);
                    break;
                }
                case '?': {
                    help();
                    break;
                }
                    return 1;
                default: {
                    abort();
                    break;
                }
            }
        }

        SDGString filename1, filename2;

        filename1 = argv[optind];
        if (++optind < argc)
            filename2 = argv[optind];

        /* Print any remaining command line arguments (not options). */
        if (++optind < argc) {
            help();
            std::cout << "non-option ARGV-elements: " << std::endl;
            while (optind < argc)
                std::cout << argv[optind++] << std::endl;
            return 1;
        }

        if (filename2.empty()) {
            repeat = true;
            filename2 = filename1;
            std::cout << "De novo mode!" <<  std::endl;
        }
        connect_dist=(kmer_dist + 1) * kmer_size;

        show_parameter(filename1, filename2);

        //Check parameters
        if (count_cutoff < 0 || count_cutoff > 1) {
            std::cout << "count_cutoff must be in interval [0-1]! Entered value is " << count_cutoff << std::endl;
            exit(EXIT_FAILURE);
        }
        if (diversity_cutoff < 0 || diversity_cutoff > 1) {
            std::cout << "diversity_cutoff must be in interval [0-1]! Entered value is " << diversity_cutoff
                      << std::endl;
            exit(EXIT_FAILURE);
        }
        if (kmask < 2 ) {
            std::cout << "kmask must be greater than 1! Entered value is " << kmask
                      << std::endl;
            exit(EXIT_FAILURE);
        }

        if (stat_only) {
            std::cout << "\nCompute kmer stat only!" << std::endl;
            for (unsigned bw = 1; bw <= bkmer_size; bw++) {
                HashDNASeq hsrch(kmer_size, kmask, mask_hole_length, algorithm,bw, kmer_dist, 0, min_size, step_q);
                std::vector<unsigned> kmer_count((unsigned) pow(4, hsrch.getEffectiveKmerSize()), 0);
                std::list<Info_kmer> list_infokmer;
                Info_kmer kmer_threshold;
                unsigned nb_kmer;

                std::cout << "\n======Compute kmer background probability for " << bw - 1
                          << " Markov's chain order======" << std::endl;

                hsrch.kmer_analysis(filename2, kmer_size, kmask, mask_hole_length, bw, kmer_size / 2, count_cutoff, diversity_cutoff,
                                    kmer_count, nb_kmer, list_infokmer, kmer_threshold);
            }
            std::cout << "\nEnd HashDNASeq (version " << VERSION << ")" << std::endl;
            exit(EXIT_SUCCESS);
        }

        Hasher hsrch(kmer_size, kmask, mask_hole_length, bkmer_size, kmer_dist, 0, min_size, step_q, pen_join, algorithm);
        //bool valid_idx_file = true;
        bool valid_idx_file = false; // suppress kidx file
        double prev_genome_perc_coverage = 0.0;
        std::stringstream alignout_name, seqout_name ;
        std::list< RangePair > frag_list;
        for (unsigned iter = 1; iter <= nb_iter || nb_iter == 0; iter++) {
            bool first_iter=false;
            if (iter==1){
                first_iter=true;
            }
            hsrch.load(filename2, kmer_size, kmask, mask_hole_length, bkmer_size, kmer_size / 2, count_cutoff, diversity_cutoff,
                       min_count,valid_idx_file, first_iter);


            FastaIstream in(filename1);
            if (!in) {
                std::cerr << "file:" << filename1 << " does not exist!" << std::endl;
                return 1;
            }

            unsigned genome_size=0;
            unsigned genome_coverage=0;
            min_frag_size = min_size;
            unsigned numseq = 0;
            std::list< RangePair > rev_frag_list, compfrag_list, rev_compfrag_list;
            frag_list.clear();

            while (in) {
                BioSeq seq;
                if (in)
                    in >> seq;
                numseq++;
                std::cout << seq.header << " len:" << seq.length() << " read!" << std::endl;
                BioSeq seq_rev=seq.reverse();
                BioSeq seq_comp=seq.complement();
                BioSeq seq_revcomp=seq_comp.reverse();
                genome_size+=seq.length();

                //Search on chunked input sequence
                if (chunk_size_kb != 0) {
                    unsigned start = 0;
                    unsigned chunk_size = ((chunk_size_kb * 1000) / kmer_size) * kmer_size;
                    unsigned nb_chunk = seq.length() / chunk_size;
                    for (unsigned i = 1; i < nb_chunk; i++) {
                        std::cout << "==>chunk #" << i << "/" << nb_chunk << ":" << start << ".."
                                  << start + chunk_size - 1 << std::endl;
                        search_frag(hsrch,seq,seq_rev,seq_comp,seq_revcomp,start,start + chunk_size - 1,numseq,
                                    connect_dist,min_frag_size,repeat,
                                    frag_list,rev_frag_list,compfrag_list,rev_compfrag_list,verbosity);
                        start = start + chunk_size;
                    }
                    std::cout << "==>chunk #" << nb_chunk << "/" << nb_chunk << ":" << start << ".." << seq.length()
                              << std::endl;
                    search_frag(hsrch,seq,seq_rev,seq_comp,seq_revcomp,start,seq.length()-1,numseq,
                                connect_dist,min_frag_size,repeat,
                                frag_list,rev_frag_list,compfrag_list,rev_compfrag_list,verbosity);
                } else {
                    //Search on full input sequence
                    search_frag(hsrch,seq,seq_rev,seq_comp,seq_revcomp,0,seq.length()-1,numseq,
                                connect_dist,min_frag_size,repeat,
                                frag_list,rev_frag_list,compfrag_list,rev_compfrag_list,verbosity);
                }
                translate_comp(frag_list,compfrag_list,seq.length());
                compfrag_list.clear();
                translate_comp(rev_frag_list,rev_compfrag_list,seq.length());
                rev_compfrag_list.clear();
                std::cout << "ok!\n" << std::endl;
            }
            in.close();


            //Compute stat on resuts and filters false positives
            unsigned qval_score,qval_len;
            double qtile=1-p_value;

            std::cout<<"--Random fragment stats:"<<std::endl;
            std::cout << "--Compute score and identity" << std::endl;


            qval_len=Hasher::fragLengthStat(rev_frag_list, qtile);
            Hasher::fragLenFilter(rev_frag_list,qval_len);
            Hasher::fragSeqAlign(rev_frag_list,filename1,filename2,true,verbosity);

            qval_score=Hasher::fragScoreStat(rev_frag_list, qtile, genome_coverage);
            std::cout<<"Coverage="<<genome_coverage<<" ("<<(float)genome_coverage/genome_size<<")"
                     <<" coverage % difference="<<fabs(((float)genome_coverage/genome_size)-prev_genome_perc_coverage)<<std::endl;

            std::cout<<"--Real fragment stats:"<<std::endl;
            std::cout << "--Compute score and identity" << std::endl;

            Hasher::fragLenFilter(frag_list,qval_len);
            Hasher::fragSeqAlign(frag_list,filename1,filename2,false,verbosity);

            Hasher::fragScoreFilter(frag_list,qval_score);
            Hasher::fragScoreStat(frag_list, qtile, genome_coverage);

            if(pen_join>0.0){
                std::cout << "Join fragment with penality " << pen_join << " ..." << std::flush;
                hsrch.fragJoin(frag_list);
                std::cout<<" done!"<<std::endl;
            }
            genome_coverage=Hasher::fragCoverage(frag_list);
            std::cout<<"**Coverage="<<genome_coverage<<" ("<<(float)genome_coverage/genome_size<<")"
                     <<" coverage % difference="<<fabs(((float)genome_coverage/genome_size)-prev_genome_perc_coverage)<<std::endl;

            //Write results
            std::ofstream alignout;
            FastaOstream seqout;
            alignout_name.str("");
            alignout_name.clear();
            seqout_name.str("");
            seqout_name.clear();
            if (!outfilename.empty()){
                alignout_name << outfilename << "." << iter << ".align";
                seqout_name << outfilename << "." << iter << ".fa";
            }
            else{
                alignout_name << filename1 << "." << iter << ".hasher.align";
                seqout_name << filename1 << "." << iter << ".hasher.fa";
            }

            std::cout<<"--Write coord in "<<alignout_name.str()<<" ... "<<std::flush;
            alignout.open(alignout_name.str());
            Hasher::fragAlignWrite(frag_list, filename1, filename2, alignout);
            alignout.close();
            std::cout<<"done!"<<std::endl;

            std::cout<<"--Write fasta in "<<seqout_name.str()<<" ... "<<std::flush;
            seqout.open(seqout_name.str());
            Hasher::fragMergeSeqWrite(frag_list, filename1, seqout);
            seqout.close();
            std::cout<<"done!"<<std::endl;

            if(genome_coverage==0) break;

            if((fabs(((float)genome_coverage/genome_size)-prev_genome_perc_coverage)<0.01 || (float)genome_coverage/genome_size>=1.0)
                && nb_iter==0 && iter>1)
                break;
            prev_genome_perc_coverage=(float)genome_coverage/genome_size;
            filename2=seqout_name.str();

            min_count++;
        }
        std::cout << "\nEnd Hasher (version " << VERSION << ")" << std::endl;
        end = clock();
        time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
        std::cout << "====>Total time spent****: " << time_spent << std::endl;

        //Write final results
        std::stringstream alignout_final_name;
        std::stringstream seqout_final_name;
        if (!outfilename.empty()) {
            alignout_final_name << outfilename << ".final.align";
            seqout_final_name << outfilename << ".final.fa";
        } else {
            alignout_final_name << filename1 << ".final.hasher.align";
            seqout_final_name << filename1 << ".final.hasher.fa";
        }

        std::stringstream cmd1,cmd2;

        cmd1<<"cp "<<alignout_name.str()<<" "<<alignout_final_name.str();
        std::system(cmd1.str().c_str());

        FastaOstream seqout;
        std::cout<<"--Write final fasta in "<<seqout_final_name.str()<<" ... "<<std::flush;
        seqout.open(seqout_name.str());
        Hasher::fragSeqWrite(frag_list, filename1, seqout);

        exit(EXIT_SUCCESS);
    }
    catch (const SDGException& e) {
        std::cerr << "******Exception catched: " << e.message << " ******" << std::endl;
        exit(EXIT_FAILURE);
    }
    catch (const std::exception &e) {
        std::cout << "Caught exception \"" << e.what() << "\"\n";
    }
    catch (const char* msg) {
        std::cerr << msg << std::endl;
        exit(EXIT_FAILURE);
    }
    catch (...) {
        std::cerr << "****** unknown exception catch !!! ******" << std::endl;
        exit(EXIT_FAILURE);
    }
}












