/*************************************************************************
 * Author: Manuel Dieterle, 2012
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMFunctor                                                            //
//                                                                      //
// Class containing all functions.                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef myTMFunctor
#define myTMFunctor

#include "TFile.h"
#include "TKey.h"
#include "TGaxis.h"
#include "TPad.h"
#include "TLatex.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TMObjectCollection.h"

namespace TMFunctor
{
   extern Bool_t fIsFunctorInit;
   extern Bool_t fReject;
   extern Bool_t fPar;
   extern Double_t fRejLo;
   extern Double_t fRejUp;
   extern Double_t parOne;
   extern Double_t parTwo;
   extern Int_t fNumber;
   extern TH1** hFit;
   extern TH2** hFit2D;

   void InitFunctor(Int_t n, TH1** h, Bool_t isReject = kFALSE, Double_t fLo = 0., Double_t fUp = 0.);
   void InitFunctor(Int_t n, TH2** h, Bool_t isReject = kFALSE, Double_t fLo = 0., Double_t fUp = 0.);
   void SetParameters(Double_t* p);
   Double_t MyFitFuncMC(Double_t *x, Double_t *par);
   Double_t MyFitFuncMC2D(Double_t *x, Double_t *par);
   Double_t FitTH1Pol1(Double_t *x, Double_t *par);
   Double_t Pol1Exc(Double_t *x, Double_t *par);
   Double_t MyFitFuncEmpty(Double_t *x, Double_t *par);
}

#endif

