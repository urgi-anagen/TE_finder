#include <FastaIstream.h>
#include <BioSeq.h>
#include <string>
#include <sstream>

FastaIstream &FastaIstream::operator>>(BioSeq &p) {
    std::ostringstream buffstr;
    bool found = false;
    std::string titre;
    std::string line;
    if (this->is_open()) {
        while (*this) {
            if (peek() == '>') {
                if (!found) {
                    if (std::getline(*this, line)) {
                        line.erase(0,1);
                        line.erase(line.find_last_not_of("\t\n\v\f\r ") + 1);
                        titre=line;
                        buffstr.str("");
                        buffstr.clear();
                    }
                    found = true;
                } else
                    break;
            } else {
                if (std::getline(*this, line))
                    buffstr << line;
            }
        }
    }
    p = buffstr.str();
    p.header=titre;
    return *this;
}

