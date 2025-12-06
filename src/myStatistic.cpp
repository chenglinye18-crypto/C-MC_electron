#include "myStatistic.h"

myStatistic::myStatistic()
{
  reset();
}

void myStatistic::reset() {
	//ntotbh
    int iband,jband;
	
    for(iband=0;iband<MNB;iband++)
    {
        NumTetChange[iband]=0;
        NumBZChange[iband]=0;
        NumCellHit[iband]=0;

        ntotp[iband]=0;
        nreap[iband]=0;
        nslfp[iband]=0;

        nreaii[iband]=0;
        dtotp[iband]=0.0;
    }

    for(iband=0;iband<NBE;iband++)
        for(jband=0;jband<MNScaEle*NBE;jband++)
            nsctype[jband][iband]=0;

    for(iband=0;iband<NBH;iband++)
        for(jband=0;jband<MNScaHole*NBH;jband++)
            nsctyph[jband][iband]=0;

    for(iband=0;iband<NBOE;iband++)
        for(jband=0;jband<MNScaOxEle*NBOE;jband++)
            nsctypoe[jband][iband]=0;
}

void Contact::reset() {
  NumParGen = 0;
  CharGen = 0;    
  NumParCatch = 0;
  CharCatch = 0;
  EnergyGen = 0;
  EnergyCatch = 0;

  NumEleGen = NumHoleGen = 0;
  NumEleCatch = NumHoleCatch = 0;
  CharEleGen = CharHoleGen = 0.0;
  CharEleCatch = CharHoleCatch = 0.0;
}
