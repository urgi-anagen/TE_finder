#ifndef DUSTER_H
#define DUSTER_H

#include <FastaIstream.h>
#include <FastaOstream.h>
#include <fstream>
#include <string>
#include <cmath>
#include <list>
#include <tuple>
#include <unordered_map>
#include "HashDNASeq.h"


class Test_Hasher;


class Duster : public HashDNASeq {
    friend class Test_Duster;


public:

    explicit Duster(unsigned w = 10, unsigned msk = 100, unsigned mask_hole_length = 1, unsigned bw = 2,
                    unsigned wd = 1, unsigned fd = 1, unsigned minsize = 20, unsigned step = 1) :
                    HashDNASeq(w, msk, mask_hole_length, 1, 1, bw, wd, fd, minsize, step) {};

    unsigned getEffectiveKmerSize() override { return hseq.getEffectiveKmerSize(); };

    void search(const BioSeq &sequence, unsigned numseq, unsigned start, unsigned end, bool repeat,
                std::vector<std::pair<unsigned, unsigned> > &fmerged, unsigned verbose);

    static void
    writeBED(const SDGString &qname, const std::vector<std::pair<unsigned, unsigned> > &frag, std::ostream &out);

    static unsigned compute_coverage(const std::vector<std::pair<unsigned, unsigned> > &frag);

    void fragMerge(std::vector<std::pair<unsigned, unsigned> > &frag,
                   unsigned connect_dist,
                   std::vector<std::pair<unsigned, unsigned> > &fmerged);

    static void get_sequences(const std::vector<std::pair<unsigned, unsigned> > &frag, BioSeq &seq, FastaOstream &out);
};


#endif






