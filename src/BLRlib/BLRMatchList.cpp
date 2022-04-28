/**
 *
 * BLRMatchList.cpp
 *
 */

#include <BLRMatchList.h>

BLRMatchList::~BLRMatchList() {
//! -- destructor --
    clear();
}

void BLRMatchList::clear() {
//! -- set empty match_list  --
    match_list.clear();
}

void BLRMatchList::read_blast_results(const std::string &match_filename, double max_eval, unsigned min_length,
                                      double min_identity, bool is_wublast) {
    std::ifstream fin(match_filename);           /*file where input matchList*/
    if (fin.bad()) {
        std::cerr << "ERROR: " << match_filename << " could not be open!" << std::endl;
        exit(EXIT_FAILURE);
    }
    bool empty = true;

    while (fin) {
        if (is_wublast) {
            try {
                match.readwublastfield(fin);
            }
            catch (BlastMatch::Empty_Parser_output e) {
                if (empty)
                    std::cout << "WuBlast empty output" << std::endl;
                break;
            }
        } else {
            try {
                match.readncbiblastfield(fin);
            }
            catch (BlastMatch::Empty_Parser_output e) {
                if (empty)
                    std::cout << "NCBIBlast empty output" << std::endl;
                break;
            }
        }
        empty = false;
        if (match.getE_value() < max_eval
            && match.getIdentity() > min_identity
            && match.getLength() > min_length) {
            match_list.push_back(match);
        }
    }
    fin.close();
}

void BLRMatchList::remove_self_hits(unsigned cut_length, unsigned cut_over) {
    //  -------------------------------------------------
    //! - Function to remove_self_hits the matches that do not need
    /*  - Arguments: length is length of cut-out
        -            over is the overlap of cut-out
        -            same_bank_name is boolen, true if bank and querybank are same
        ---------------------------------------------------*/

    std::list<BlastMatch>::iterator iter_list = match_list.begin();  /*iterator to build list of matches*/

    while (iter_list != match_list.end()) {
        bool remove = false;
        // Removing alignement of the same query/subject
        if ((iter_list->getQuery_num() == iter_list->getSubject_num())
            && (iter_list->getQuery_strt() == iter_list->getSubject_strt())
            && (iter_list->getQuery_end() == iter_list->getSubject_end())) {
            remove = true;
        } else if (
                ((iter_list->getQuery_num() == iter_list->getSubject_num() - 1)
                 &&
                 (iter_list->getQuery_strt() == (cut_length - cut_over + 1)
                  && iter_list->getQuery_end() == cut_length)
                 &&
                 (iter_list->getSubject_strt() == 1
                  && iter_list->getSubject_end() == cut_length)
                 &&
                 (cut_length > 0 && cut_over > 0))
                ||
                ((iter_list->getQuery_num() - 1 == iter_list->getSubject_num())
                 &&
                 (iter_list->getSubject_strt() == (cut_length - cut_over + 1)
                  && iter_list->getSubject_end() == cut_length)
                 &&
                 (iter_list->getQuery_strt() == 1
                  && iter_list->getQuery_end() == cut_over)
                 &&
                 (cut_length > 0 && cut_over > 0))
                ) {
            remove = true;
        }
        if (remove)
            iter_list = match_list.erase(iter_list);
        else
            iter_list++;
    }
}

void BLRMatchList::show_list() {
    // --------------------------------------------------
    //!- Function for developement to show the matchlist
    /* --------------------------------------------------*/
    std::list<BlastMatch>::iterator id = match_list.begin();
    std::cout << "-----------------" << std::endl;
    while (id != match_list.end()) {
        id->view();
        id++;
    }
}

void BLRMatchList::save_list(SDGString name, int flag) {
    // -----------------------------------------------------
    /*!- Function to save the matchlist
      /param  name is file to save
      /param flag is boolean true if is a new run of blast*/
    //-----------------------------------------------------
    std::ofstream file;
    if (flag)
        // to write in file end
        file.open(name.start(), std::ios::app);
    else
        // to write in file begin
        file.open(name.start(), std::ios::trunc);
    if (file.bad()) {
        std::cout << "ERROR: " << name.start() << " could not open" << std::endl;
        return;
    }

    savebin(file);

    file.close();
}

//! fonction to save the hash_align in binary file
void BLRMatchList::savebin(std::ostream &out) {
    std::list<BlastMatch>::iterator iter_list = match_list.begin();  //iterator to visit match_list
    while (iter_list != match_list.end()) {
        iter_list->writelst(out);
        iter_list++;
    }
}

//! fonction to write the hash_align in ascii format
void BLRMatchList::writetxt(std::ostream &out) {
    std::list<BlastMatch>::iterator iter_list = match_list.begin();  //iterator to visit match_list
    while (iter_list != match_list.end()) {
        out << *iter_list;
        iter_list++;
    }
}

//! fonction to read the hash_align in ascii format
void BLRMatchList::readtxt(std::istream &in) {
    BlastMatch match;            // current alignment
    while (in) {
        in >> match;
        if (match.getQuery_num() == 0) break;
        match_list.push_back(match);
    }
}

void BLRMatchList::readtxt(char *data) {
    std::string sdata = data;
    std::istringstream sin(sdata);
    BlastMatch match;            // current alignment
    while (sin) {
        sin >> match;
        if (match.getQuery_num() == 0) break;
        match_list.push_back(match);
    }
}
