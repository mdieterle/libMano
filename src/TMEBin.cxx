/*************************************************************************
 * Author: Manuel Dieterle, 2012
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMEBin                                                               //
//                                                                      //
// Class representing an energy bin.                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TMEBin.h"

ClassImp(TMEBin)

//______________________________________________________________________________
TMEBin::TMEBin(Int_t nTC, Int_t iTC, Int_t TClo, Int_t TCup, Double_t TCe, Double_t TCee)
    : TNamed()
{
    // Constructor. Allocate space for 'nMaxTaggCh' tagger channels.
 
    // init members
    fEnergy = TCe;
    fEnergyError = TCee;
    nTaggCh = nTC;
    iTaggCh = iTC;
    fNTaggCh = 1;
    fTaggChLo = TClo;
    fTaggChUp = TCup;

    for (Int_t i = fTaggChLo; i < fTaggChUp; i++) fNTaggCh++;
}

//______________________________________________________________________________
TMEBin::~TMEBin()
{
    // Destructor.
}

//______________________________________________________________________________
void TMEBin::Print()
{
    // Print out the content of this class.

    printf("TMEBin content:\n");
    printf("Energy                : %f\n", fEnergy);
    printf("Energy error          : %f\n", fEnergyError);
    printf("total # of tagg ch.   : %d\n", nTaggCh);
    printf("number of tagg ch.    : %d\n", iTaggCh);
    printf("# of tagger ch.       : %d\n", fNTaggCh);
    printf("first tagger ch.      : %d\n", fTaggChLo);
    printf("last tagger ch.       : %d\n", fTaggChUp);
}

