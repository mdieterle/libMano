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


#ifndef myTMEBin
#define myTMEBin

#include "TNamed.h"

class TMEBin : public TNamed
{

protected:
    Double_t fEnergy;                           // energy
    Double_t fEnergyError;                      // energy error
    Int_t nTaggCh;                              // total number of tagger channels
    Int_t iTaggCh;                              // number of tagger channel
    Int_t fNTaggCh;                             // number of tagger channels contributing to this bin
    Int_t fTaggChLo;                            // lowest tagger channel contributing to this bin
    Int_t fTaggChUp;                            // highest tagger channel contributing to this bin

public:
    TMEBin() : TNamed(), fEnergy(0.), fEnergyError(0.), nTaggCh(0),
               iTaggCh(0), fNTaggCh(0), fTaggChLo(0), fTaggChUp(0) { }
    TMEBin(Int_t nTaggCh, Int_t iTaggCh, Int_t TClo, Int_t TCup, Double_t TCe, Double_t TCee);
    virtual ~TMEBin();
    void Print();

    Double_t GetEnergy() const { return fEnergy; }
    Double_t GetEnergyError() const { return fEnergyError; }
    Int_t GetNTotalTaggerChannels() const { return nTaggCh; }
    Int_t GetNTaggerChannel() const { return iTaggCh; }
    Int_t GetNTaggerChannels() const { return fNTaggCh; }
    Int_t GetTaggerChannelLo() const { return fTaggChLo; }
    Int_t GetTaggerChannelUp() const { return fTaggChUp; }
    
    ClassDef(TMEBin, 1)  // Bin of differential cross section
};

#endif

