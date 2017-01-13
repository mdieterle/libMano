/*************************************************************************
 * Author: Manuel Dieterle, 2012
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMEnergyThetaCut                                                     //
//                                                                      //
// Mano Addon for TOEnergyThetaCut from OSCAR                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TMEnergyThetaCut.h"

ClassImp(TMEnergyThetaCut)

//______________________________________________________________________________
TMEnergyThetaCut::TMEnergyThetaCut(Int_t nCTBins, const Char_t* name)
    : TNamed(name, name)
//    : TOEnergyThetaCut(nCTBins, name)
{
    // Constructor using 'nCTBins' cos(theta) bins.
    // The object will be named 'name'.
    // set the members
    fNCTBin = nCTBins;

    // init arrays
    fCutLowerFunc = new TF1*[fNCTBin];
    fCutUpperFunc = new TF1*[fNCTBin];

    for (Int_t i = 0; i < fNCTBin; i++)
    {
        fCutLowerFunc[i] = 0;
        fCutUpperFunc[i] = 0;
    }
}

//______________________________________________________________________________
TMEnergyThetaCut::~TMEnergyThetaCut()
{
    if (fCutLowerFunc)
    {
        for (Int_t i = 0; i < fNCTBin; i++)
        {
            if (fCutLowerFunc[i]) delete fCutLowerFunc[i];
            if (fCutUpperFunc[i]) delete fCutUpperFunc[i];
        }
        delete [] fCutLowerFunc;
        delete [] fCutUpperFunc;
    }
}

//______________________________________________________________________________
void TMEnergyThetaCut::SetLowerFunc(Int_t ctBin, TF1* f)
{
    // Set the lower function for the cos(theta) bin 'ctBin'. The function 'f' will
    // be copied into this class.

    Char_t tmp[256];

    // check index range
    if (ctBin < 0 || ctBin >= fNCTBin)
    {
        Error("SetLowerFunc", "cos(theta) bin %d cannot be found! (number of cos(theta) bins: %d)", ctBin, fNCTBin);
        return;
    }

    // clone graph
    sprintf(tmp, "Lower_%d", ctBin);
    fCutLowerFunc[ctBin] = (TF1*) f->Clone(tmp);
    fCutLowerFunc[ctBin]->SetTitle(tmp);
}

//______________________________________________________________________________
void TMEnergyThetaCut::SetUpperFunc(Int_t ctBin, TF1* f)
{
    // Set the upper function for the cos(theta) bin 'ctBin'. The function 'f' will
    // be copied into this class.

    Char_t tmp[256];

    // check index range
    if (ctBin < 0 || ctBin >= fNCTBin)
    {
        Error("SetUpperFunc", "cos(theta) bin %d cannot be found! (number of cos(theta) bins: %d)", ctBin, fNCTBin);
        return;
    }

    // clone graph
    sprintf(tmp, "Upper_%d", ctBin);
    fCutUpperFunc[ctBin] = (TF1*) f->Clone(tmp);
    fCutUpperFunc[ctBin]->SetTitle(tmp);
}

//______________________________________________________________________________
Bool_t TMEnergyThetaCut::IsInsideS(Double_t energy, Double_t cosTheta, Double_t value)
{
    // Return kTRUE if the value 'value' for the energy 'energy' and cos(theta)
    // value 'cosTheta' lies within the functions.

    // check cos(theta) value
    if (cosTheta < -1 || cosTheta >= 1)
    {
        Error("IsInsideS", "Bad cos(theta) parameter value! (%f)", cosTheta);
        return kFALSE;
    }

    // calculate the cos(theta) bin
    Int_t ctBin = Int_t(fNCTBin * (cosTheta + 1.) / 2.);

    // get mean and sigma
    Double_t lower, upper;
    if (fCutLowerFunc[0])
    {
        lower = fCutLowerFunc[ctBin]->Eval(energy);
        upper = fCutUpperFunc[ctBin]->Eval(energy);
    }
    else
    {
       Error("IsInsideS","No lower and upper functions found");
       return kFALSE;
    }

    if (value > lower && 
        value < upper) return kTRUE;
    else return kFALSE;
}
//______________________________________________________________________________
