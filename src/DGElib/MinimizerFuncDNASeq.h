//
// Created by Hadi Quesneville on 03/01/2022.
//

#ifndef TE_FINDER_MINIMIZERFUNCDNASEQ_H
#define TE_FINDER_MINIMIZERFUNCDNASEQ_H

struct MinimizerFuncDNASeq // Compute hashing value of a word
{
    unsigned kmer_size;
    unsigned window_size;



    MinimizerFuncDNASeq(unsigned k, unsigned w) : kmer_size(k), window_size(w) {
        if (kmer_size > 15)
            throw SDGException(NULL, "MinimizerFuncDNASeq: kmer size must be < 16 !!");
    };

    virtual ~MinimizerFuncDNASeq() {
    }
    unsigned getEffectiveKmerSize(void){return kmer_size;};
    unsigned minimizer(const char *p, unsigned pos_start, unsigned &pos) {

        const char *w = p;
        unsigned val, min_val = 0, min_h = 0;
        for (unsigned i = 0; i < window_size && *w != '\0'; i++, w++) {
            unsigned h = hash(w, val);
            if (i == 0 || val < min_val) {
                min_h = h;
                min_val = val;
                pos=pos_start+i;
            } else if(val == min_val && h < min_h)
            {
                min_h = h;
                pos=pos_start+i;
            }
        }
        return min_h;
    }

    unsigned hash(const char *p, unsigned & val) {
        unsigned h = 0, val_nuc = 0;
        val=0;
        for (unsigned i = 0; i < kmer_size && *p != '\0'; i++, p++) {
            switch (*p) {

                case 'C': {
                    val_nuc = 3;
                    val+=1;
                    break;
                }
                case 'G': {
                    val_nuc = 2;
                    val+=1;
                    break;
                }
                case 'T': {
                    val_nuc = 1;
                    break;
                }
                case 'c': {
                    val_nuc = 3;
                    val+=1;
                    break;
                }
                case 'g': {
                    val_nuc = 2;
                    val+=1;
                    break;
                }
                case 't': {
                    val_nuc = 1;
                    break;
                }
                default: {
                    val_nuc = 0;
                    break;
                }
            }
            h <<= 2;
            h |= val_nuc;
        }
        return h;
    };

    unsigned operator()(const char *p, unsigned pos_start, unsigned &pos) { return minimizer(p,pos_start,pos); };
};
#endif //TE_FINDER_MINIMIZERFUNCDNASEQ_H
