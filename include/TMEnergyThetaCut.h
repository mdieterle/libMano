
/*************************************************************************
 * Author: Manuel Dieterle, 2012
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMEnergyThetaCut                                                         //
//                                                                      //
// Class for a collection of TObjects                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef myTMEnergyThetaCut
#define myTMEnergyThetaCut

#include "TF1.h"
//#include "TOEnergyThetaCut.h"

//class TMEnergyThetaCut : public TOEnergyThetaCut
class TMEnergyThetaCut : public TNamed
{

private:
    Int_t fNCTBin;                              // number of cos(theta) bins
    TF1** fCutLowerFunc;                        //[fNCTBin] cut lower function
    TF1** fCutUpperFunc;                        //[fNCTBin] cut upper function

public:
//    TMEnergyThetaCut() : TOEnergyThetaCut(),
    TMEnergyThetaCut() : TNamed(),
                         fNCTBin(0),
                         fCutLowerFunc(0), fCutUpperFunc(0) { }
    TMEnergyThetaCut(Int_t nCTBins, const Char_t* name);
    virtual ~TMEnergyThetaCut();

    void SetLowerFunc(Int_t ctBin, TF1* f);
    void SetUpperFunc(Int_t ctBin, TF1* f);
    TF1* GetLowerFunction(Int_t ctBin) const { return fCutLowerFunc[ctBin]; }
    TF1* GetUpperFunction(Int_t ctBin) const { return fCutUpperFunc[ctBin]; }

    Bool_t IsInsideS(Double_t energy, Double_t cosTheta, Double_t value);

    ClassDef(TMEnergyThetaCut, 1)  // Bin of differential cross section
};

#endif

