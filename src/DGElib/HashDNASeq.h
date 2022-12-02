#ifndef HASHDNASEQ_H
#define HASHDNASEQ_H

#include <SDGError.h>
#include <SDGString.h>
#include <BioSeq.h>
#include <FastaIstream.h>
#include <fstream>
#include <string>
#include <cmath>
#include <list>
#include <set>
#include <tuple>
#include <unordered_map>
#include <map>
#include <type_traits>

#include "HashFuncDNASeq.h"
#include "MinimizerFuncDNASeq.h"


class Test_HashDNASeq;


struct Info_kmer {
    unsigned count;
    double expected_count;
    double entropy;
    double diversity;
    double good_kmer;
    unsigned hash_key;

    explicit Info_kmer(unsigned h = 0, unsigned c = 0, double ec = 0, double e = 0, double d = 0, double g = 0) :
            count(c),
            expected_count(ec),
            entropy(e),
            diversity(d),
            good_kmer(g),
            hash_key(h) {};

    bool operator<(const Info_kmer &k) const {
        if (count < k.count) return true;
        return false;
    };
};

struct Compare_Info_kmer {
    bool operator()(const Info_kmer &k1, const Info_kmer &k2) {
        return k1 < k2;
    };
};


class HashDNASeq {
    friend class Test_HashDNASeq;

    friend class Test_Hasher;

protected:

    struct KmerSpos //Store a kmer position in a sequence
    {
        unsigned pos;
        unsigned numSeq;

        explicit KmerSpos(unsigned p = 0, unsigned n = 0) : pos(p), numSeq(n) {
          if (n>100000000)
              throw Unsigned_Out_Of_Range("KmerSpos - sequence number setup out of range: n=",n);
          if(p>1000000000)
              throw Unsigned_Out_Of_Range("KmerSpos - sequence position setup out of range: p=",p);
            pos = p;
            numSeq = n;
        };

        friend int operator<(const KmerSpos &w1, const KmerSpos &w2) {
            if (w1.numSeq < w2.numSeq) return 1;
            if (w1.numSeq > w2.numSeq) return 0;
            if (w1.numSeq == w2.numSeq && w1.pos < w2.pos) return 1;
            return 0;
        };

        friend int operator==(const KmerSpos &w1, const KmerSpos &w2) {
            if (w1.numSeq != w2.numSeq) return 0;
            if (w1.pos != w2.pos) return 0;
            return 1;
        };

        friend std::ostream &operator<<(std::ostream &os, const KmerSpos &w) {
            os << "(" << w.numSeq << "," << w.pos << ")";
            return os;
        };

    };

    struct Diag // Store the diagonal of a kmer match
    {
        long diag;
        KmerSpos wpos;

        explicit Diag(long d = 0, unsigned p = 0, unsigned n = 0) {
            diag = d;
            wpos = KmerSpos(p, n);
        };

        friend std::ostream &operator<<(std::ostream &os, const Diag &d) {
            os << "[" << d.wpos << "," << d.diag << "]";
            return os;
        };

        friend int operator<(const Diag &d1, const Diag &d2) {
            if (d1.wpos.numSeq < d2.wpos.numSeq) return 1;
            if (d1.wpos.numSeq > d2.wpos.numSeq) return 0;
            if (d1.wpos.numSeq == d2.wpos.numSeq) {
                if (d1.diag < d2.diag) return 1;
                if (d1.diag == d2.diag)
                    if (d1.wpos.pos < d2.wpos.pos) return 1;
            }
            return 0;
        };
    };

/*    struct Diag_map : public std::vector< std::unordered_map<long,std::list<Diag> > >{

            Diag_map(unsigned size){this->resize(size+1);};
            void insert(unsigned numseq, long diag, Diag val){
                this->at(numseq)[diag].emplace_back(val);
            }
    };*/
    typedef std::vector< std::list<Diag> > Diag_map;

    HashFuncDNASeq hseq;
    MinimizerFuncDNASeq mseq;
    HashFuncDNASeq bhseq, mhseq, nhseq;

    std::unordered_map<unsigned, std::vector<KmerSpos>::iterator> hash2wpos, hash_ptr;
    std::vector<KmerSpos> kmer_pos;

    std::vector<SDGString> subjectName;
    unsigned kmer_size, window_size, bkmer_size, mkmer_size, wdist, fdist, min_size, step_q, max_key, nfrag;
    unsigned nbseqQ, nbseqS;
    unsigned hash_algorithm;

    struct Key : std::pair<long, long> {
        Key(long i, long j) {
            first = i;
            second = j;
        };
    };

    void kmer_counts(const SDGString &filenameS, unsigned kmerSize, unsigned bkmerSize, unsigned mkmerSize,
                     std::vector<unsigned> &wcount, unsigned &nb_kmer,
                     std::vector<unsigned> &bcount, unsigned &nb_bkmer,
                     std::vector<unsigned> &mcount, unsigned &nb_mkmer,
                     std::vector<unsigned> &ncount, unsigned &nb_nuc);

    void kmer_prob(unsigned wsize, unsigned bwsize, unsigned mwsize, unsigned mask_hole_period, unsigned mask_hole_length,unsigned kmer_window,
                   const std::vector<unsigned> &wcount, unsigned nb_kmer,
                   const std::vector<unsigned> &bcount, unsigned nb_bkmer,
                   const std::vector<unsigned> &mcount, unsigned nb_mkmer,
                   const std::vector<unsigned> &ncount, unsigned nb_nuc,
                   std::list<Info_kmer> &list_infokmer) const;

    static void kmer_occurrence_percentiles(const std::list<Info_kmer> &list_infokmer, double cutoff,
                                            Info_kmer &kmer_threshold);

    static void kmer_entropy_percentiles(const std::list<Info_kmer> &list_infokmer, double cutoff_entropy,
                                  Info_kmer &kmer_threshold);

    static void kmer_diversity_percentiles(const std::list<Info_kmer> &list_infokmer, double cutoff_diversity,
                                    Info_kmer &kmer_threshold);

    static void kmer_goodkmer_percentiles(const std::list<Info_kmer> &list_infokmer);

    void kmer_filter(const std::list<Info_kmer> &list_infokmer, const Info_kmer &kmer_threshold, unsigned min_count,
                     std::vector<unsigned> &wcount, bool first_iter) const;

    void kmer_ssr_filter(unsigned wsize, std::vector<unsigned> &wcount);

    bool read_idx(const SDGString &filename, double count_cutoff, double diversity_cutoff, unsigned min_count,
                  unsigned mask_hole_period, unsigned mask_hole_length);

    void save_idx(const SDGString &filename, double count_cutoff, double diversity_cutoff, unsigned min_count,
                  unsigned mask_hole_period, unsigned mask_hole_length, const std::vector<unsigned> &wcount);

    unsigned hashSeqCountWHole(const BioSeq &seq, unsigned wsize, std::vector<unsigned> &wcount);
    unsigned hashSeqCount(const BioSeq &seq, unsigned wsize, std::vector<unsigned> &wcount);

    void hashSubjectSeqPosWHole(const BioSeq &seq, const std::vector<unsigned> &wcount);
    void hashSubjectSeqPos(const BioSeq &seq, unsigned num_seq, unsigned wsize, const std::vector<unsigned> &wcount);
    void hashSubjectSeqPosMinimizer(const BioSeq &seq, unsigned num_seq, unsigned wsize, const std::vector<unsigned> &wcount);
    void hashSubjectSeqPosWHoleMinimizer(const BioSeq &seq, unsigned num_seq, unsigned wsize, const std::vector<unsigned> &wcount);

    void hashSeqPos(const BioSeq &seq, unsigned wsize, std::list<std::pair<unsigned,unsigned>>& kmer_pos_list);
    void hashSeqPosWHole(const BioSeq &seq, unsigned wsize, std::list<std::pair<unsigned,unsigned>>& kmer_pos_list);

    void matchKmersHole(const BioSeq &sequence, unsigned qSeq, unsigned start, unsigned end, bool repeat,
                        Diag_map &diag_map);
    void matchKmers(const BioSeq &sequence, unsigned qSeq, unsigned start, unsigned end, bool repeat,
                    Diag_map &diag_map);
    void matchKmersMinimizer(const BioSeq &sequence, unsigned qSeq, unsigned start, unsigned end, bool repeat,
                             Diag_map &diag_map);
    void matchKmersWHoleMinimizer(const BioSeq &sequence, unsigned qSeq, unsigned start, unsigned end, bool repeat,
                             Diag_map &diag_map);

    void minimize(unsigned window_size, std::list<std::pair<unsigned, unsigned>> kmer_pos_list,
                  std::set<std::pair<unsigned, unsigned>> &minimized_kmer_pos_list);

    static void diagSearchDist(Diag_map& diag_map,
                        unsigned connect_dist, unsigned kmerSize,
                        std::vector<std::pair<unsigned, unsigned> > &frag);

public:

    HashDNASeq(unsigned kmer_size = 10, unsigned mask_hole_period = 10, unsigned mask_hole_length = 1,
               unsigned kmer_window=20, unsigned alg=1, unsigned bw = 2, unsigned wd = 1,
               unsigned fd = 1, unsigned minsize = 20, unsigned step = 1) :
            hseq(kmer_size, mask_hole_period, mask_hole_length), mseq(kmer_size, kmer_window),
            bhseq(bw, 0, 0), mhseq(kmer_size / 2, 0, 0), nhseq(1, 0, 0),
            kmer_size(kmer_size), window_size(kmer_window),bkmer_size(bw), mkmer_size(kmer_size / 2),
            wdist(wd),
            fdist(fd),
            min_size(minsize),
            step_q(step),
            nbseqQ(0),
            nbseqS(0),
            hash_algorithm(alg){
        if(alg==2){
            hseq.set(kmer_size,0,0);
        }
        if (hseq.getEffectiveKmerSize() > 16)
            throw SDGException(NULL, "HashDNASeq: Kmer size must be <= 16 !!");
        max_key = (unsigned) pow(4, hseq.getEffectiveKmerSize());
    };

    virtual unsigned getEffectiveKmerSize() { return hseq.getEffectiveKmerSize(); };

    virtual void
    load(const SDGString &filenameS, unsigned kmerSize, unsigned mask_hole_period, unsigned mask_hole_length, unsigned kmer_window, unsigned bkmerSize,
         unsigned mkmerSize, double count_cutoff, double diversity_cutoff,
         unsigned min_count,
         bool &valid_idx_file, bool first_iter, bool filter_ssr);

    void kmer_analysis(const SDGString &filenameS,
                       unsigned kmerSize, unsigned mask_hole_period, unsigned mask_hole_length,
                       unsigned kmer_window,
                       unsigned bkmerSize, unsigned mkmerSize,
                       double count_cutoff, double diversity_cutoff,
                       std::vector<unsigned> &wcount,unsigned &nb_kmer,
                       std::list<Info_kmer> &list_infokmer, Info_kmer &kmer_threshold);

    void search(const BioSeq &sequence, unsigned numseq, unsigned start, unsigned end, bool repeat,
                std::vector<std::pair<unsigned, unsigned> > &frag, unsigned verbose=0);
};


#endif






