/*
 * \file Cutter.ccp
 */

#include <limits.h>
#include <Cutter.h>
#include <SDGMemBioSeq.h>
#include <SDGFastaIstream.h>
#include <SDGFastaOstream.h>
#include <RangeSeq.h>
#include <RangeMap.h>
#include <SDGBioSeqDB.h>

bool operator==(const Cutter &op1, const Cutter &op2) {
    if (op1.length == op2.length &&
        op1.over == op2.over &&
        op1.word == op2.word &&
        op1.extention == op2.extention)
        return true;
    else
        return false;
}

bool operator!=(const Cutter &op1, const Cutter &op2) {
    if (op1.length != op2.length ||
        op1.over != op2.over ||
        op1.word != op2.word ||
        op1.extention != op2.extention)
        return true;
    else
        return false;
}

SDGString Cutter::cutDB(SDGString bank_name, int verbose) {
    SDGFastaIstream cutin(bank_name);
    if (!cutin) {
        std::cerr << "ERROR: wrong file name '" << bank_name << "'";
        return extention;
    }
    if (verbose > 0)
        std::cout << "Cutting bank '" << bank_name << "':" << std::endl;

    SDGString cutfile;
    cutfile = bank_name + extention;

    if (length < over) {
        std::cerr << "ERROR: cut-out length must be greater than overlap" << std::endl;
        return extention;
    }

    SDGFastaOstream cutout(cutfile);
    if (!cutout) {
        std::cerr << "ERROR: could not open '" << cutfile << "'";
        return extention;
    }

    RangeMap clrseq_map;
    RangeMap range_mapN;
    RangeSeq range;
    if (verbose > 0) {
        if (word == 0)
            std::cout << "parsing sequences..." << std::endl;
        else
            std::cout << "parsing sequences and reading low informative regions..." << std::endl;
    }

    SDGFastaIstream bank_in(bank_name);
    if (!bank_in) {
        std::cerr << "file:" << bank_name << " does not exist!" << std::endl;
    }
    while (bank_in) {
        SDGBioSeq seq;
        if (!bank_in)  break;
        bank_in >> seq;
        if (verbose > 0) {
            std::cout << seq.getDE() << ": " << seq.length() << " bp" << std::flush;
            if (seq.length() == 0)
                std::cout << " --> Cannot cut empty sequence, skip it!" << std::endl;
            else
                std::cout << " " <<seq.toString() << std::endl;
        }

        if (seq.length() > 0)
                clrseq_map.add(RangeSeq("", seq.getDE(), 1, seq.length()));

        unsigned long start = 0, pos = 0;
        unsigned nbN = 0;
        unsigned last_pos = 0;
        while (pos < seq.length() && word != 0) {
            while (pos < seq.length()
                   && (seq.charAt(pos) == 'N' || seq.charAt(pos) == 'n'
                       || seq.charAt(pos) == 'X' || seq.charAt(pos) == 'x')) {
                if (nbN == 0)
                    start = pos;
                nbN++;
                pos++;
                last_pos = pos;
            }
            if (pos - last_pos > word) {
                if (nbN > word) {
                    if (verbose > 0)
                        std::cout << "Found N_stretch "<<seq.getDE()<<":"<<(start + 1)<<".."<<last_pos << std::endl;
                    range.set("N_stretch", seq.getDE(), start + 1, last_pos);
                    range_mapN.add(range);
                }
                nbN = 0;
            }
            pos++;
        }
        if (nbN > word) {
            if (verbose > 0)
                std::cout << "Found N_stretch "<<seq.getDE()<<":"<<(start + 1)<<".."<<last_pos << std::endl;
            range.set("N_stretch", seq.getDE(), start + 1, last_pos);
            range_mapN.add(range);
        }
    }
    bank_in.close();

    range_mapN.save(bank_name + ".Nstretch.map");

    if (word != 0) {
        if (verbose > 0)
            std::cout << "subtracting low informative region..." << std::endl;
        clrseq_map.diff(range_mapN);
    }

    if (length > 0) {
        if (verbose > 0)
            std::cout << "cutting long fragments..." << std::endl;
        clrseq_map.cut(length, over);
    }

    if (verbose > 0)
        std::cout << "writing cut bank..." << std::endl;
    clrseq_map.writeCutSeq(cutfile, bank_name, verbose - 1);

    if (verbose > 0)
        std::cout << "Bank '" << bank_name << "' was cut." << std::endl;

    return cutfile;
}

bool Cutter::check(SDGString bank_name, int verbose) {
    SDGString cutfile = bank_name + extention;
    if (verbose > 0)
        std::cout << "Checking bank '" << cutfile << "'...";

    SDGFastaIstream cutin(cutfile);
    if (!cutin) {
        if (verbose > 0)
            std::cout << " not cut !" << std::endl;
        return false;
    }

    if (length < over) {
        std::cerr << "ERROR: cut length less than overlap" << std::endl;
        return false;
    }
    if (length == 0) {
        std::cout << "WARNING: cut length is zero" << std::endl;
    }
    if (length > ULONG_MAX) {
        std::cerr << "ERROR: cut length is too big" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (verbose > 0)
        std::cout << "file '" << cutfile << "' already exists" << std::endl;
    return true;
}
