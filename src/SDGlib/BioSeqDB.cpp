//
// Created by Hadi Quesneville on 29/04/2022.
//

#include "BioSeqDB.h"
#include "FastaIstream.h"
#include "SDGError.h"
#include <fstream>


//-------------------------------------------------------
void BioSeqDB::load(std::string fichier, int verbose) {
    clear();
    filename = fichier;
    if (verbose > 0)
        std::cout << "Loading file '" << fichier << "'... " << std::endl;
    FastaIstream file(fichier);
    unsigned count = 0;

    while (file) {
        BioSeq s;
        file >> s;
        push_back(s);
        if (name2pos.find(back().header) != name2pos.end()) {
            std::ostringstream ostr;
            ostr << "SDGBioSeqDB error: Duplicated entry '" << back().header << "' in " << fichier;
            throw SDGException(NULL, ostr.str(), -1);
        }
        name2pos[back().header] = count++;
    }
    if (verbose > 0)
        std::cout << "Files was loaded." << std::endl;
}
