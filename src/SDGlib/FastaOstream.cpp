#include <FastaOstream.h>




FastaOstream &FastaOstream::operator << (BioSeq p)
{
    std::ofstream &out = *this;
    out << ">" << p.header << std::endl;

    for (unsigned long l = 0; l < p.size(); l += 60)
        out << p.substr(l, 60) << "\n";

    return *this;
}














