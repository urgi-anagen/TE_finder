/**
 * \file RangeAlign.h
 * \brief Header file for the class RangeAlign
 */

#ifndef RANGEALIGN_H
#define RANGEALIGN_H

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <string>
#include <SDGError.h>
#include "Range.h"

/**
 * \class RangeAlign
 * \brief Implement a fragment
 *
 * This class inherits from Range and adds a name and number to the fragment.
 */
class RangeAlign: public Range
{
  long num_chr;
  std::string name_seq;

	protected:

	public:
		RangeAlign(const RangeAlign *o);
		RangeAlign(const RangeAlign& o);
		RangeAlign(long c=0, ulong s=0, ulong e=0);
		RangeAlign(std::string name_seq,long c=0, ulong s=0, ulong e=0);
		RangeAlign(const Range& r)
		{
			num_chr=-1;
			name_seq="";
			Range::set(r.getStart(),r.getEnd());
		};

		virtual ~RangeAlign(void);
		virtual void *clone(void) const;


		void set(long c=0, ulong s=0, ulong e=0)
		{
			name_seq="";
			num_chr=c;
			Range::set(s,e);
		};

		void set(std::string name,long c=0, ulong s=0, ulong e=0)
		{
			name_seq=name;
			num_chr=c;
			Range::set(s,e);
		};

        void translate_comp(unsigned len_seq){
            Range::translate_comp(len_seq);
        }

		long getNumChr(){return num_chr;};
		long getNumChr() const {return num_chr;};
		void setNumChr(long c){num_chr=c;};

		void setNameSeq(std::string name) {name_seq=name;};
		std::string getNameSeq(void) {return name_seq;};
		std::string getNameSeq(void) const {return name_seq;};

		RangeAlign diff(const RangeAlign& r)
		{
			Range rd=Range::diff(r);
			return RangeAlign(r.name_seq,r.num_chr,rd.getStart(),rd.getEnd());
		};

		bool overlap(const RangeAlign& r) const
		{
			if(num_chr!=r.num_chr || name_seq!=r.name_seq) return false;
			return Range::overlap(r);
		};

		bool overlap(const RangeAlign& r)
		{
			if(num_chr!=r.num_chr || name_seq!=r.name_seq) return false;
			return Range::overlap(r);
		};

		long distance(const RangeAlign& r)
		{
			if(num_chr!=r.num_chr || name_seq!=r.name_seq)
				throw SDGException(this,
						"RangeAlign::distance: Not on same chromosome!");
			return Range::distance(r);
		};

		long distance(const RangeAlign& r) const
		{
			if(num_chr!=r.num_chr || name_seq!=r.name_seq)
				throw SDGException(this,
						"RangeAlign::distance: Not on same chromosome!");
			return Range::distance(r);
		};

		void view(void) const
		{
		std::cout<<"num chr: "<<num_chr<<" name seq: "<<name_seq<<std::endl;
		Range::view();
		}

		friend bool operator<(const RangeAlign& r1, const RangeAlign& r2)
		{
			if(std::abs(r1.num_chr)<std::abs(r2.num_chr)
					|| r1.name_seq<r2.name_seq) return true;
			else if(std::abs(r1.num_chr)==std::abs(r2.num_chr) &&  r1.name_seq==r2.name_seq)
				return operator<((Range)r1,(Range)r2);
			return false;
		};

		friend bool operator==(const RangeAlign& r1,const RangeAlign& r2)
		{
			if(r1.num_chr!=r2.num_chr || r1.name_seq!=r2.name_seq) return false;
			else return operator==((Range)r1,(Range)r2);
		};

		friend bool operator<=(const RangeAlign& r1, const RangeAlign& r2)
		{
			if(r1<r2 || r1==r2) return true;
			return false;
		};

		friend std::ostream& operator<<(std::ostream& out,const RangeAlign& r)
		{
			out<<r.name_seq<<'\t'<<r.num_chr<<"\t"<<r.start<<"\t"<<r.end;
			return out;
		};

		friend std::istream& operator>>(std::istream& in, RangeAlign& r)
		{
			in>>r.name_seq>>r.num_chr>>r.start>>r.end;
			return in;
		};
};
const RangeAlign emptyRangeAlign;

#endif /* RANGEALIGN_H */
