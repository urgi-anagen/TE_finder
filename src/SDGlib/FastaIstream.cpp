#include <FastaIstream.h>
#include <BioSeq.h>
#include <string>
#include <sstream>

FastaIstream &FastaIstream::operator>>(BioSeq &p) {
    std::ostringstream buffstr;
    std::string titre;
    std::string line;
    if (this->is_open()) {
        while(std::getline(*this, line)) {
            if (line[0] == '>') {
                line.erase(0, 1);
                line.erase(line.find_last_not_of("\t\n\v\f\r ") + 1);
                p.clear();
                p.header = line;
            } else {
                p+=line;
            }
            if (peek() == '>' || peek() == EOF) {
                break;
            }
        }
    }
    return *this;
}

