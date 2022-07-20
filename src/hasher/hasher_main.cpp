#include <getopt.h>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <list>

#include <chrono>
using namespace std::chrono;

#include "SDGString.h"
#include "FastaOstream.h"

#include "Hasher.h"


double filter_cutoff=0.0, pen_join=0.5;
unsigned kmer_size=15, kmer_window=20, step_q=15, bkmer_size=1, kmer_dist=5,mask_hole_length=1,connect_dist=1,
min_frag_size=20, chunk_size_kb=0, min_count=0, mask_hole_period=4, verbosity=0, overlap=0, nb_iter=1,
algorithm=1,align_ext=0, min_identity=65, min_score=25;

double count_cutoff=1.0, diversity_cutoff=0.0;
bool repeat=false, stat_only=false;
double p_value=0.05;

SDGString outfilename="";

void help(void) {
    std::cerr << "usage: hasher"
              << " [<options>] <fasta query sequence> [<fasta sequence model>]." << std::endl
              << " options:" << std::endl
              << "   -h, --help:\n\t this help" << std::endl
              << "   -w, --kmer:\n\t kmer length, default: " << kmer_size << std::endl
              << "   -W, --kmer_window:\n\t kmer window for minimizer: " << kmer_window << std::endl
              << "   -S, --step_q:\n\t step on query sequence, default: " << step_q << std::endl
              << "   -k, --mask_hole_period:\n\t period of k-mer hole, default: " << mask_hole_period << std::endl
              << "   -l, --len_hole_mask:\n\t kmer mask hole length, default: " << mask_hole_length << std::endl
              << "   -d, --kmer_dist:\n\t max number of kmer between two matching kmer to connect, default: "
              << kmer_dist << std::endl
              << "   -j, --pen_join:\n\t penality to join two matching kmer, default: "
              << pen_join << " (enter a value < 0.0 to skip this)" << std::endl
              << "   -s, --min_frag_size:\n\t min frag size to store, default: "
              << min_frag_size << std::endl
              << "   -C, --filter_cutoff:\n\t filter kmer with counts over a percentile (Value [0-1]), default: "
              << count_cutoff << std::endl
              << "   -D, --diversity_cutoff:\n\t filter kmer with diversity measure of kmer size used for background probability (Value [0-1]), default: "
              << diversity_cutoff << std::endl
              << "   -m, --min_count:\n\t filter kmer with counts less than this value, default: "
              << min_count << std::endl
              << "   -b, --background_kmer_size:\n\t kmer size to compute kmer background probability, default: "
              << bkmer_size << std::endl
              << "   -i, --min_identity:\n\t minimum fragment identity , default: " << min_identity << std::endl
              << "   -x, --min_score:\n\t minimum fragment score , default: " << min_score << std::endl
              << "   -p, --p-value:\n\t p-value, default: " << p_value << std::endl
              << "   -o, --file_out:\n\t filename for output," << std::endl
              << "   -c, --chunk_size:\n\t sequence chunk size in kb, default: None" << std::endl
              << "   -n, --nb_iter:\t number of iteration: " << nb_iter << std::endl
              << "   -a, --analysis:\n\t compute kmer statistics only" << std::endl
              << "   -A, --algorithm:\n\t algorithm number (0: simple kmers, 1: kmers with hole, 2: kmer minimizer, 3: kmers with hole and minimizer), default:"
              << algorithm << std::endl
              << "   -e, --extend:\n\t extend initial hit by alignment (0: no extension, 1: full SW alignment, 2: align on boundaries ), default:"
              << align_ext << std::endl
              << "   -v, --verbosity:\n\t verbosity level, default:" << verbosity << std::endl;
};
void show_parameter(SDGString filename1,SDGString filename2) {
    std::cout << "\nRun with parameters:\n"
              << "Query sequences: " << filename1 << std::endl
              << "Model sequences: " << filename2 << std::endl
              << "   -w, --kmer_size:\t kmer length: " << kmer_size << std::endl
              << "   -W, --kmer_window:\t kmer window for minimizer: " << kmer_window << std::endl
              << "   -S, --step_q:\t step on query sequence: " << step_q << std::endl
              << "   -k, --mask_hole_period:\t period of k-mer hole: " << mask_hole_period << std::endl
              << "   -l, --len_hole_mask:\t kmer mask hole length: " << mask_hole_length << std::endl
              << "   -d, --kmer_dist:\t max number of kmer between two matching kmer to connect: " << kmer_dist
              << std::endl
              << "   -j, --pen_join:\t penality to join two matching kmer, default: "
              << pen_join << " (enter a value < 0.0 to skip this)" << std::endl
              << "   -s, --min_frag_size:\t min frag size to store: " << min_frag_size << std::endl
              << "   -C, --filter_cutoff:\t filter kmer with counts in the last percentile: " << count_cutoff
              << std::endl
              << "   -D, --diversity_cutoff:\t filter kmer with diversity measure of kmer size used for background probability: "
              << diversity_cutoff << std::endl
              << "   -m, --min_count:\t filter kmer with counts less than this value: " << min_count << std::endl
              << "   -b, --background_kmer_size:\t kmer size to compute kmer background probability: " << bkmer_size
              << std::endl
              << "   -i, --min_identity:\t minimum fragment identity , default: " << min_identity << std::endl
              << "   -x, --min_score:\t minimum fragment score , default: " << min_score << std::endl
              << "   -p, --p-value:\t p-value, default: " << p_value << std::endl
              << "   -o, --file_out:\t filename for output:" << outfilename << std::endl
              << "   -c, --chunk_size:\t sequence chunk size in kb: " << chunk_size_kb << std::endl
              << "   -n, --nb_iter:\t number of iteration: " << nb_iter << std::endl
              << "   -A, --algorithm:\t algorithm number (0: simple kmers, 1: kmers with hole, 2: kmer minimizer, 3: kmers with hole and minimizer ) : " << algorithm
              << std::endl
              << "   -e, --extend:\t extend initial hit by alignment (0: no extension, 1: full SW alignment, 2: align on boundaries ) :"
              << align_ext << std::endl
              << "   -v, --verbosity:\t verbosity level: " << verbosity << std::endl;
};
// search on sequence chunk, reverse sequence, and reverse complement
void search_frag(Hasher& hsrch,
                 const BioSeq& seq, const BioSeq& seq_rev,
                 const BioSeq& seq_comp, const BioSeq& seq_revcomp,
                 unsigned start, unsigned end, unsigned numseq, unsigned Kmer_connect_dist,
                 unsigned min_fragment_size, bool denovo_mode,
                 std::list< RangePair >& frag_list,
                 std::list< RangePair >& rev_frag_list,
                 std::list< RangePair >& compfrag_list,
                 std::list< RangePair >& rev_compfrag_list,
                 double p_value,
                 unsigned verbosity)
{
    if(verbosity>0)
        std::cout<<"++> run search on direct sequence."<<std::endl;
    hsrch.search(seq, start, end, numseq, Kmer_connect_dist,
                 min_fragment_size, denovo_mode, frag_list, verbosity);
    if(p_value != 1.0){
        if(verbosity>0)
            std::cout<<"++> run search on direct reversed (random) sequence."<<std::endl;
        hsrch.search(seq_rev, start, end, numseq, Kmer_connect_dist,
                     min_fragment_size, denovo_mode, rev_frag_list, verbosity);
    }

    unsigned rev_start=seq_comp.size()-1-end;
    unsigned rev_end=seq_comp.size()-1-start;
    if(verbosity>0)
        std::cout<<"++> run search on reverse complementary sequence."<<std::endl;
    hsrch.search(seq_comp, rev_start, rev_end, numseq, Kmer_connect_dist,
                 min_fragment_size, denovo_mode, compfrag_list, verbosity);
    if(p_value != 1.0){
        if(verbosity>0)
            std::cout<<"++> run search on reverse complementary reversed (random) sequence."<<std::endl;
        hsrch.search(seq_revcomp, rev_start, rev_end, numseq, Kmer_connect_dist,
                     min_fragment_size, denovo_mode, rev_compfrag_list, verbosity);
    }
}
//translate reverse complement coordinate
void translate_comp(std::list< RangePair >& frag_list,const std::list< RangePair >& compfrag_list,
                    unsigned len_seq){
    for (auto rp : compfrag_list){
        rp.getRangeQ().translate_comp(len_seq);
        frag_list.push_back(rp);
    }
}

//substract from query, ranges that were already found
void range_substract(unsigned num_chr, Range query_range, const std::list<RangePair>& match_list, std::list<Range>& substracted_query_range) {
    substracted_query_range.emplace_back(query_range);
    for(std::list<RangePair>::const_iterator it_m=match_list.begin(); it_m!=match_list.end(); it_m++){
        if(it_m->getRangeQ().getNumChr()==num_chr){
            for(std::list<Range>::iterator it_q=substracted_query_range.begin(); it_q!=substracted_query_range.end(); it_q++){
                Range new_range=it_q->diff(it_m->getRangeQ());
                if(!new_range.empty()) substracted_query_range.emplace_back(new_range);
                if(it_q->empty() || it_q->getLength()<min_frag_size) it_q=substracted_query_range.erase(it_q);
            }
        }
    }
}

int main(int argc, char *argv[]) {
    try {
        std::cout << "Beginning Hasher (version " << VERSION << ")" << std::endl;
        auto begin = high_resolution_clock::now();

        if (argc == 1) {
            help();
            exit(EXIT_SUCCESS);
        }
        int c;
        while (1) {
            static struct option long_options[] =
                    {
                            {"help",                  no_argument,       0, 'h'},
                            {"kmer",                  required_argument, 0, 'w'},
                            {"kmer_window",           required_argument, 0, 'W'},
                            {"step_q",                required_argument, 0, 'S'},
                            {"mask_hole_period",                 required_argument, 0, 'k'},
                            {"len_hole_mask",         required_argument, 0, 'l'},
                            {"kmer_dist",             required_argument, 0, 'd'},
                            {"pen_join",              required_argument, 0, 'j'},
                            {"min_frag_size",              required_argument, 0, 's'},
                            {"filter_cutoff",         required_argument, 0, 'C'},
                            {"min_count",            required_argument, 0, 'm'},
                            {"diversity_cutoff",     required_argument, 0, 'D'},
                            {"background_kmer_size", required_argument, 0, 'b'},
                            {"min_identity",              required_argument, 0, 'i'},
                            {"min_score",              required_argument, 0, 'x'},
                            {"p_value",              required_argument, 0, 'p'},
                            {"file_out",             required_argument, 0, 'o'},
                            {"chunk_size",           required_argument, 0, 'c'},
                            {"nb_iter",              required_argument, 0, 'n'},
                            {"analysis",             no_argument,       0, 'a'},
                            {"algorithm",            required_argument, 0, 'A'},
                            {"extend",            required_argument, 0, 'e'},
                            {"verbosity",            no_argument,       0, 'v'},
                            {0, 0,                                      0, 0}
                    };
            /* `getopt_long' stores the option index here. */
            int option_index = 0;

            c = getopt_long(argc, argv, "hd:f:w:W:S:k:l:d:j:s:C:D:m:b:i:x:p:o:c:n:aA:e:v:",
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
                case 'W': {
                    kmer_window = atoi(optarg);
                    break;
                }
                case 'S': {
                    step_q = atoi(optarg);
                    break;
                }
                case 'k': {
                    mask_hole_period = atoi(optarg);
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
                    min_frag_size = atoi(optarg);
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
                case 'i': {
                    min_identity = atoi(optarg);
                    break;
                }
                case 'x': {
                    min_score = atoi(optarg);
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
                case 'e': {
                    align_ext = atoi(optarg);
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
        if (mask_hole_period < 2 ) {
            std::cout << "mask_hole_period must be greater than 1! Entered value is " << mask_hole_period
                      << std::endl;
            exit(EXIT_FAILURE);
        }

        if (stat_only) {
            std::cout << "\nCompute kmer stat only!" << std::endl;
            for (unsigned bw = 1; bw <= bkmer_size; bw++) {
                HashDNASeq hsrch(kmer_size, mask_hole_period, mask_hole_length, algorithm, bw, kmer_dist, 0, min_frag_size, step_q);
                std::vector<unsigned> kmer_count((unsigned) pow(4, hsrch.getEffectiveKmerSize()), 0);
                std::list<Info_kmer> list_infokmer;
                Info_kmer kmer_threshold;
                unsigned nb_kmer;

                std::cout << "\n======Compute kmer background probability for " << bw - 1
                          << " Markov's chain order======" << std::endl;

                hsrch.kmer_analysis(filename2, kmer_size, mask_hole_period, mask_hole_length, kmer_window,
                                    bkmer_size,(kmer_size / 2), count_cutoff, diversity_cutoff,
                                    kmer_count, nb_kmer, list_infokmer, kmer_threshold);
            }
            std::cout << "\nEnd HashDNASeq (version " << VERSION << ")" << std::endl;
            exit(EXIT_SUCCESS);
        }

        //bool valid_idx_file = true;
        bool valid_idx_file = false; // suppress kidx file
        double prev_genome_perc_coverage = 0.0;
        std::stringstream alignout_name, seqout_name ;
        std::list< RangePair > frag_list, all_frag_list;
        for (unsigned iter = 1; iter <= nb_iter || nb_iter == 0; iter++) {
            bool first_iter=false;
            if (iter==1){
                first_iter=true;
            }
            else if (iter==2){
                min_count++;
            }

            Hasher hsrch(kmer_size, mask_hole_period, mask_hole_length, kmer_window, bkmer_size, kmer_dist, 0, min_frag_size, step_q, pen_join, algorithm);
            if(hsrch.getEffectiveKmerSize()<7){
                std::cout << "\nkmer effective size = " << hsrch.getEffectiveKmerSize() << " < 7, stopping search !" << std::endl;
                break;
            }
            hsrch.load(filename2, kmer_size, mask_hole_period, mask_hole_length, kmer_window, bkmer_size, kmer_size / 2, count_cutoff, diversity_cutoff,
                       min_count, valid_idx_file, first_iter);


            FastaIstream in(filename1);
            if (!in) {
                std::cerr << "file:" << filename1 << " does not exist!" << std::endl;
                return 1;
            }

            unsigned genome_size=0;
            unsigned genome_coverage=0;
            unsigned query_length=0;
            unsigned numseq = 0;
            std::list< RangePair >  rev_frag_list, compfrag_list, rev_compfrag_list;
            frag_list.clear();

            while (in) {
                //Load query sequence
                BioSeq seq;
                if (in)
                    in >> seq;
                numseq++;
                std::cout << seq.header << " len:" << seq.length() << " read!" << std::endl;
                BioSeq seq_rev=seq.reverse();
                BioSeq seq_comp=seq.complement();
                BioSeq seq_revcomp=seq_comp.reverse();
                genome_size+=seq.length();

                Range query_range(1,seq.length());
                std::list<Range> query_range_list;
                range_substract(numseq, query_range, all_frag_list, query_range_list);


                unsigned chunk_size = ((chunk_size_kb * 1000) / kmer_size) * kmer_size;
                for(std::list<Range>::const_iterator it_r=query_range_list.begin();it_r!=query_range_list.end();it_r++){
                    std::cout << "=>Query #" << numseq << ": " << it_r->getStart() << ".." << it_r->getEnd() << std::endl;
                    unsigned start = it_r->getStart()-1;
                    unsigned end = it_r->getEnd()-1;
                    query_length+=(end-start+1);

                    //Search on chunked input sequence
                    if (chunk_size_kb != 0) {
                        unsigned nb_chunk = it_r->getLength() / chunk_size;
                        for (unsigned i = 1; i < nb_chunk; i++) {
                            if(verbosity>0)
                                std::cout << "==>chunk #" << i << "/" << nb_chunk << ":" << start + 1 << ".."
                                      << start + chunk_size << std::endl;
                            search_frag(hsrch,seq,seq_rev,seq_comp,seq_revcomp,start,start + chunk_size - 1,numseq,
                                        connect_dist,min_frag_size,repeat,
                                        frag_list,rev_frag_list,compfrag_list,rev_compfrag_list,p_value,verbosity);
                            start = start + chunk_size;
                        }
                        if(end-start>min_frag_size){
                            if(verbosity>0)
                                std::cout << "==>chunk #" << nb_chunk << "/" << nb_chunk << ":" << start + 1
                                << ".." << end + 1 << std::endl;
                            search_frag(hsrch,seq,seq_rev,seq_comp,seq_revcomp,start,end,numseq,
                                        connect_dist,min_frag_size,repeat,
                                        frag_list,rev_frag_list,compfrag_list,rev_compfrag_list,p_value,verbosity);
                        }
                    } else {
                        //Search on full input sequence
                        search_frag(hsrch,seq,seq_rev,seq_comp,seq_revcomp,start,start + it_r->getLength()-1,numseq,
                                    connect_dist,min_frag_size,repeat,
                                    frag_list,rev_frag_list,compfrag_list,rev_compfrag_list,p_value,verbosity);
                    }
                }

                translate_comp(frag_list,compfrag_list,seq.length());
                compfrag_list.clear();
                if(p_value!=1.0){
                    translate_comp(rev_frag_list,rev_compfrag_list,seq.length());
                    rev_compfrag_list.clear();
                }
                std::cout << "ok!\n" << std::endl;
            }
            in.close();
            unsigned qval_score=min_score, qval_identity = min_identity, qval_len=0;
            double qtile=0.95;
            if(p_value!=1.0) {
                //Compute stat on resuts and filters false positives
                qtile = 1 - p_value;

                std::cout << "--Random fragment stats: iteration "<<iter<<" / kmer size:"<<kmer_size<<std::endl;
                std::cout << "--Compute score and identity" << std::endl;


                qval_len = Hasher::fragLengthStat(rev_frag_list, qtile);
                Hasher::fragLenFilter(rev_frag_list, qval_len);
                Hasher::fragSeqAlign(rev_frag_list, filename1, filename2, true, verbosity);

                qval_score = Hasher::fragScoreIdentityStat(rev_frag_list, qtile, genome_coverage);
                unsigned rev_genome_coverage= Hasher::fragCoverage(rev_frag_list);
                std::cout << "Coverage=" << genome_coverage << " (" << (float) genome_coverage / query_length << ")"
                          << " merged coverage % = "
                          << (float) rev_genome_coverage / query_length
                          << std::endl;
            }
            std::cout<<"--Real fragment stats: iteration "<<iter<<" / kmer size:"<<kmer_size<<std::endl;
            std::cout << "--Compute score and identity" << std::endl;

            Hasher::fragLenFilter(frag_list,qval_len);
            Hasher::fragSeqAlign(frag_list,filename1,filename2,false,verbosity);

            Hasher::fragScoreIdentityFilter(frag_list, qval_score, qval_identity);
            Hasher::fragScoreIdentityStat(frag_list, qtile, genome_coverage);

            if(pen_join>0.0){
                std::cout << "Join fragment with penality " << pen_join << " ..." << std::flush;
                hsrch.fragJoin(frag_list);
                std::cout<<" done!"<<std::endl;
            }

            std::cout << "Align fragments" << " ..." <<std::endl;
            unsigned ext_len=20;
            if(align_ext==1){
                hsrch.fragSeqSWAlign(frag_list, ext_len, filename1, filename2, false, verbosity);
            }
            if(align_ext==2){
                hsrch.fragSeqExtAlign(frag_list, ext_len, filename1, filename2, false, verbosity);
            }

            std::cout<<" done!"<<std::endl;

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

            //concat coords
            std::list<RangePair> copy_frag_list=frag_list;
            all_frag_list.splice(all_frag_list.end(),copy_frag_list);
            genome_coverage=Hasher::fragCoverage(all_frag_list);
            std::cout<<"**Coverage iteration "<<iter<<" / kmer size:"<<kmer_size<<" = "<<genome_coverage<<" ("<<(float)genome_coverage/genome_size<<")"
                     <<" coverage % difference="<<fabs(((float)genome_coverage/genome_size)-prev_genome_perc_coverage)<<std::endl;

            if(genome_coverage==0) break;

            if((fabs(((float)genome_coverage/genome_size)-prev_genome_perc_coverage)<0.01 || (float)genome_coverage/genome_size>=1.0)
                && nb_iter==0 && iter>1)
                break;
            prev_genome_perc_coverage=(float)genome_coverage/genome_size;

            //concat align files
            if(iter==1) remove(alignout_final_name.str().c_str());
            std::ofstream oalign_a(alignout_final_name.str(), std::ios_base::binary | std::ios_base::app);
            std::ifstream ialign_b(alignout_name.str(), std::ios_base::binary);
            oalign_a.seekp(0, std::ios_base::end);
            oalign_a << ialign_b.rdbuf();

            //concat fasta files
            if(iter==1) remove(seqout_final_name.str().c_str());
            std::ofstream oseq_a(seqout_final_name.str(), std::ios_base::binary | std::ios_base::app);
            std::ifstream iseq_b(seqout_name.str(), std::ios_base::binary);
            oseq_a.seekp(0, std::ios_base::end);
            oseq_a << iseq_b.rdbuf();

            filename2=seqout_final_name.str();
            if(kmer_size>1) kmer_size--;
            if(mask_hole_period>1) mask_hole_period--;
            if(step_q>1) step_q--;
            auto finish = high_resolution_clock::now();
            duration<double> elapsed = finish - begin;
            std::cout << "==>Elapsed time at iteration "<<iter<<" : " << elapsed.count() / 60 << " min." << std::endl;
        }
        std::cout << "\nEnd Hasher (version " << VERSION << ")" << std::endl;
        auto finish = high_resolution_clock::now();
        duration<double> elapsed = finish - begin;
        std::cout << "====>Total time spent****: " << elapsed.count() / 60 << " min."<<std::endl;

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
/*    catch (...) {
        std::cerr << "****** unknown exception catch !!! ******" << std::endl;
        exit(EXIT_FAILURE);
    }*/
}













