#ifndef ITRALIGN_H
#define ITRALIGN_H

#include "ExtAlign.h"

class ItrAlign: public ExtAlign
{
  unsigned max_align_len;
  double id_thres;
  unsigned min_len;
  unsigned startseq2,endseq2;
    
 public:

  ItrAlign(unsigned len, double id, unsigned max_al_len=500) : ExtAlign()
    {
      max_align_len=max_al_len;
      id_thres=id;
      min_len=len;
    };

  bool hasItr(SDGBioSeq seq);
  void alignItr(SDGBioSeq seq);

 int getStartSeq1(){return starti;};
 int getStartSeq2(){return startseq2;};
 int getEndSeq1(){return endi;};
 int getEndSeq2(){return endseq2;};
};
#endif
