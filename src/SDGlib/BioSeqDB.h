//
// Created by Hadi Quesneville on 29/04/2022.
//

#ifndef TE_FINDER_BIOSEQDB_H
#define TE_FINDER_BIOSEQDB_H

#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include "BioSeq.h"

class BioSeqDB : public std::vector<BioSeq> {
    std::string filename;
    std::map<std::string, unsigned> name2pos;
public:

    BioSeqDB() = default;
    explicit BioSeqDB(std::string
                      fichier,
                      int verbose = 0
    ) {
        load(fichier, verbose
        );
    };

    void load(std::string fichier, int verbose = 0);

    BioSeq &find(std::string name) { return operator[](name2pos[name]); };

};


#endif //TE_FINDER_BIOSEQDB_H
