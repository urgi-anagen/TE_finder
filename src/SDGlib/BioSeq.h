//
// Created by Hadi Quesneville on 22/04/2021.
//

#ifndef TE_FINDER_BIOSEQ_H
#define TE_FINDER_BIOSEQ_H
#include <string>
#include <sstream>
#include <stdexcept>

class BioSeq : public std::string {

public:

    std::string header;

    BioSeq(const std::string& seq="", const std::string& head="") : std::string(seq){
        header=head;
    }

    BioSeq subseq(unsigned start, unsigned len){
        BioSeq sseq;
        if(start+len<=size())
            sseq=substr(start,len);
        else
            sseq=substr(start);
        std::ostringstream head;
        head<<header<<"("<<start<<","<<len<<")";
        sseq.header=head.str();
        return sseq;
    }
    BioSeq complement(void){
        const std::string DNADGalphabet = "ACGTMRSVWYHKDBNXacgtmrsvwyhkdbnx";
        const std::string DNAcomp = "TGCAKYSBWRDMHVNXtgcakysbwrdmhvnx";
        std::ostringstream cstr;
        for (std::string::const_reverse_iterator i = rbegin();
             i != rend(); i++){
            std::string::size_type rep = DNADGalphabet.find(*i);
            if(rep==std::string::npos)
            {
                std::string msg="Can't find a complement to nucleotide: ";
                msg+="a";
                msg+=" !!\n";
                throw std::invalid_argument(msg);
            }
            cstr<<DNAcomp[rep];
        }
        std::ostringstream head;
        head<<header<<"(complement)";
        BioSeq cseq(cstr.str(),head.str());
        return cseq;
    }
    BioSeq reverse(void) {
        std::ostringstream ostr;
        for (std::string::const_reverse_iterator i = rbegin();
             i != rend(); i++)
            ostr << *i;
        std::ostringstream head;
        head<<header<<"(reverse)";
        BioSeq rseq(ostr.str(),head.str());
        return rseq;
    }

};


#endif //TE_FINDER_BIOSEQ_H
