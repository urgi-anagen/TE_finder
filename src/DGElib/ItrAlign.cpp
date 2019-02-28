#include "ItrAlign.h"

//---------------------------------------------------------------------------
void ItrAlign::alignItr(SDGBioSeq seq)
{
  SDGBioSeq dirseq, revseq;
  unsigned len=seq.length();
  revseq=seq.complement();
  dirseq=seq.subseq(0,(len/2<max_align_len?len/2:max_align_len));
  revseq=revseq.subseq(0,(len/2<max_align_len?len/2:max_align_len));
  setSeq(dirseq,revseq);	  
  align();
  startseq2=len-endj+1;
  endseq2=len-startj+1;
}
//---------------------------------------------------------------------------
bool ItrAlign::hasItr(SDGBioSeq seq)
{
  alignItr(seq);
  unsigned length=getEndSeq1()-getStartSeq1()+1;

  std::cout<<std::endl;
  view();

  if(getScore()==0) return false;

  std::cout<<"\n\tQuality="<<(double)getScore()/length
      <<"\n\t("<<getStartSeq1()
      <<","<<getEndSeq1()<<")"
      <<"..("<<getStartSeq2()
      <<","<<getEndSeq2()<<")"
      <<" Length itr="<<length
      <<"\n\tLength TE="<<getEndSeq2()-getStartSeq1()+1
      <<std::endl;

  if( length>=min_len && getIdentity()>id_thres )
    {
      return true;
    }
  return false;
}
