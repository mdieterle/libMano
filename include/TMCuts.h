/*************************************************************************
 * Author: Manuel Dieterle, 2012
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMCuts                                                               //
//                                                                      //
// Class containing methods to determine analysis cuts.                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef myTMCuts
#define myTMCuts

#include <iostream>
#include <fstream>
#include <stdio.h>
#include "TFile.h"
#include "TGaxis.h"
#include "TPad.h"
#include "TLatex.h"
#include "TMObjectCollection.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TMTools.h"
#include "TF1.h"
#include "TMFunctor.h"
#include "TOEnergyThetaCut.h"
#include "TOHUtils.h"

namespace TMCuts
{
    void FitMMnew(TMObjectCollection* m, const Char_t* szName, Int_t iCT = 1, Int_t iRebinMM = 10, TOEnergyThetaCut* e = 0);
    void FitMM(TMObjectCollection* m, const Char_t* szName, Int_t iCT = 1, Int_t iRebinMM = 10, Bool_t kFit = kTRUE, Bool_t kEmpty = kFALSE, Bool_t k2Empty = kFALSE);
    void SubtractIMfromMM(Int_t iCT, TFile* fDataIn, TFile* fMCIn, const Char_t* szHin, const Char_t* szHout);
    TH2F* SubtractEmptyTarget(TFile* f, TFile* fCorr, const Char_t* szName);
    void FitCopNew(TMObjectCollection* m, const Char_t* szName, Int_t iCT = 1, Int_t iRebinMM = 10, TOEnergyThetaCut* e = 0);
    void FitCop(TMObjectCollection* m, const Char_t* szName, Int_t iCT = 1, Int_t iRebinMM = 10, Bool_t kFit = kTRUE, Double_t fLeft = 100., Double_t fRight = 260.);
    void FitDEvsE(TH2F* h, Int_t iR, Double_t fSigma = 3., Bool_t kMC = kFALSE);
    void FitDEvsTOFx(TH2F* h, Int_t iR, Double_t fSigma = 3., Bool_t kMC = kFALSE);
    void FitDEvsTOFy(TH2F* h, Int_t iR);
}

#endif

