
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
            name2pos[back().getDE()] = count++;
            if(verbose>0)
                std::cout<<back().getDE()<<" length:"<<back().length()<<" ...loaded!"<<std::endl;
        }
        file.get();
    }
    if (verbose > 0)
        std::cout << "Files was loaded." << std::endl;

    file.close();
}
