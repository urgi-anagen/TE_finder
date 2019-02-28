#include <SDGFastaOstream.h>

//SDGFastaOstream SDGFastaCout(1);
//SDGFastaOstream SDGFastaCerr(2);


SDGFastaOstream &SDGFastaOstream::operator << (SDGBioSeq p)
{
  SDGString seq= ">";
  std::ofstream &out = *this;

  SDGString BQ = p.getBQ();
  SDGString AC = p.getAC();
  SDGString ID = p.getID();
  SDGString DE = p.getDE();

  if (AC.length() && ! ID.length()) ID = AC;
  if (AC.length() && ! BQ.length()) BQ = "bq";

  if (!DE.length()) 
    {
      DE += "Sequence length = " + SDGString(p.length());
    }

  if (AC.length())
    seq += BQ + "|" + AC + "|" + ID + " " + DE + "\n";
  else
    seq +=  DE + "\n";

  out << seq;

  for (unsigned long l=0; l< p.length(); l+=60)
      out << p.toString(l,60) << "\n";

  return *this;
}














