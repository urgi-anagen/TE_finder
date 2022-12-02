
#include <SDGBioSeqDB.h>
#include <SDGFastaBioSeq.h>
#include <fstream>


//-------------------------------------------------------
void SDGBioSeqDB::load(SDGString fichier, int verbose) {
    clear();
    filename = fichier;
    if (verbose > 0)
        std::cout << "Loading file '" << fichier << "'... " << std::endl;
    std::ifstream file(fichier);
    unsigned count = 0;

    while (file) {
        if (file.peek() == '>') {
            push_back(newSDGFastaBioSeq(fichier, file.tellg()));
            if (name2pos.find(back().getDE()) != name2pos.end()) {
                std::ostringstream ostr;
                ostr << "SDGBioSeqDB error: Duplicated entry '" << back().getDE() << "' in " << fichier;
                throw SDGException(NULL, ostr.str(), -1);
            }
            name2pos[back().getDE()] = count++;
            if (verbose > 0)
                std::cout << back().getDE() << " length:" << back().length() << " ...loaded!" << std::endl;
            seqlen.push_back(back().length());
        }
        file.get();
    }
    if (verbose > 0)
        std::cout << "Files was loaded." << std::endl;

    file.close();
}
