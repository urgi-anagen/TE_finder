#ifndef RANGESEQ_H
#define RANGESEQ_H

#include <iostream>
#include <algorithm>
#include <string>
#include <string.h>
#include <list>
#include <SDGError.h>
#include "Range.h"
#include "RangeAlignSet.h"

class RangeSeq : public Range {
    std::list<Range> range_set;
    std::string name, chr;

protected:

public:
    RangeSeq(const RangeSeq *o);

    RangeSeq(RangeSeq const &o);

    RangeSeq(const std::string &n = "", const std::string &c = "", ulong s = 0, ulong e = 0);

    RangeSeq(const RangeAlignSet &r,
             const std::string &n = "", const std::string &c = "");

    virtual ~RangeSeq();

    virtual void *clone() const;

    bool empty(void) { return (name == "" && chr == "") || (start == 0 && end == 0); };

    void set(const std::string &n = "", const std::string &c = "", ulong s = 0, ulong e = 0) {
        range_set.clear();
        name = n;
        chr = c;
        Range::set(s, e);
    };

    void set(const RangeAlignSet &r,
             const std::string &n = "", const std::string &c = "") {
        name = n;
        chr = c;
        Range::set(r.getStart(), r.getEnd());
        range_set = r.getRangeSet();
    };

    const std::string &getName(void) const { return name; };

    const std::string &getName(void) { return name; };

    void setName(const std::string &n) { name = n; };

    const std::string &getChr(void) const { return chr; };

    const std::string &getChr(void) { return chr; };

    void setChr(const std::string &c) { chr = c; };

    const std::list<Range> &getRangeSet(void) const { return range_set; };

    bool overlap(const RangeSeq &r) {
        if (chr != r.chr) return false;
        return Range::overlap(r);
    };

    long distance(const RangeSeq &r) {
        if (chr != r.chr)
            throw SDGException(this,
                               "RangeSeq::distance: Not on same chromosome!");
        return Range::distance(r);
    };

    RangeSeq diff(const RangeSeq &r) {
        Range rd = Range::diff(r);
        return RangeSeq(name, chr, rd.getStart(), rd.getEnd());
    };

    friend bool operator<(const RangeSeq &r1, const RangeSeq &r2) {
        if (r1.chr < r2.chr) return true;
        else return operator<((Range) r1, (Range) r2);
    };

    friend bool operator==(const RangeSeq &r1, const RangeSeq &r2) {
        if (r1.chr != r2.chr) return false;
        else return operator==((Range) r1, (Range) r2);
    };

    friend std::istream &operator>>(std::istream &in, RangeSeq &r) {
        char buff[256];
        int i, j;

        in.getline(buff, 255, '\t');
        i = 0;
        while (buff[i++] == ' ');
        j = strlen(buff);
        while (buff[--j] == ' ' && j > i);
        buff[j + 1] = '\0';
        r.name = &buff[i - 1];


        in.getline(buff, 255, '\t');
        i = 0;
        while (buff[i++] == ' ');
        j = strlen(buff);
        while (buff[--j] == ' ' && j > i);
        buff[j + 1] = '\0';
        r.chr = &buff[i - 1];

        in.getline(buff, 255, '\t');
        r.start = atoi(buff);

        in.getline(buff, 255, '\n');
        r.end = atoi(buff);

        return in;
    };

    friend std::ostream &operator<<(std::ostream &out, const RangeSeq &r) {
        out << r.name << "\t" << r.chr << "\t" << r.start << "\t" << r.end;
        return out;
    };

    void write_rangeSet(std::ostream &out, unsigned id) {
        for (std::list<Range>::const_iterator i = range_set.begin();
             i != range_set.end(); i++)
            out << id << "\t" << name << "\t" << chr << "\t"
                << i->getStart() << "\t" << i->getEnd() << std::endl;
    };
};

#endif /* RANGESEQ_H */


