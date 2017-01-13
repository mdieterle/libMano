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


#include "TMFunctor.h"

Bool_t TMFunctor::fIsFunctorInit = kFALSE;
Bool_t TMFunctor::fReject = kFALSE;
Bool_t TMFunctor::fPar = kFALSE;
Double_t TMFunctor::fRejLo = 0.;
Double_t TMFunctor::fRejUp = 0.;
Double_t TMFunctor::parOne = 0.;
Double_t TMFunctor::parTwo = 0.;
Int_t TMFunctor::fNumber;
TH1** TMFunctor::hFit = 0;
TH2** TMFunctor::hFit2D = 0;

//______________________________________________________________________________
void TMFunctor::InitFunctor(Int_t n, TH1** h, Bool_t isReject, Double_t fLo, Double_t fUp)
{
    fNumber = n;
    hFit = new TH1*[fNumber];

    for (Int_t i = 0; i < fNumber; i++)
       hFit[i] = h[i];

    if (isReject)
    {
       fRejLo  = fLo;
       fRejUp  = fUp;
       fReject = kTRUE;
    }

    fIsFunctorInit = kTRUE;

    return;
}

//______________________________________________________________________________
void TMFunctor::InitFunctor(Int_t n, TH2** h, Bool_t isReject, Double_t fLo, Double_t fUp)
{
    fNumber = n;
    hFit2D = new TH2*[fNumber];

    for (Int_t i = 0; i < fNumber; i++)
       hFit2D[i] = h[i];

    if (isReject)
    {
       fRejLo  = fLo;
       fRejUp  = fUp;
       fReject = kTRUE;
    }

    fIsFunctorInit = kTRUE;

    return;
}

//______________________________________________________________________________
void TMFunctor::SetParameters(Double_t* p)
{
   parOne = p[0];
   parTwo = p[1];

   fPar = kTRUE;

   return;
}

//______________________________________________________________________________
Double_t TMFunctor::MyFitFuncMC(Double_t *x, Double_t *par)
{
    if (!fIsFunctorInit)
       Error("TMFunctor::MyFitFuncMC","Functor not initalized!");

    Double_t iVal[fNumber];
    Double_t fVal = 0.;

    for (Int_t i = 0; i < fNumber; i++)
    {
       iVal[i] = hFit[i]->GetBinContent(hFit[i]->FindBin(x[0]));
       fVal   += par[i]*iVal[i];
    }

    return fVal;
}

//______________________________________________________________________________
Double_t TMFunctor::MyFitFuncMC2D(Double_t *x, Double_t *par)
{
    if (!fIsFunctorInit)
       Error("TMFunctor::MyFitFuncMC2D","Functor not initalized!");

    Double_t iVal[fNumber];
    Double_t fVal = 0.;

    for (Int_t i = 0; i < fNumber; i++)
    {
       iVal[i] = hFit2D[i]->GetBinContent(hFit2D[i]->FindBin(x[0],x[1]));
       fVal   += par[i]*iVal[i];
    }

    return fVal;
}

//______________________________________________________________________________
Double_t TMFunctor::FitTH1Pol1(Double_t *x, Double_t *par)
{
    if (!fIsFunctorInit)
       Error("TMFunctor::MyFitFunc","Functor not initalized!");

    Double_t iVal[fNumber];
    Double_t fVal = 0.;

    for (Int_t i = 0; i < fNumber; i++)
    {
       iVal[i] = hFit[i]->GetBinContent(hFit[i]->FindBin(x[0]));
       fVal   += par[i]*iVal[i];
    }

    if (fPar) fVal += parOne + parTwo*x[0]; 
    else Error("TMFunctor::FitTH1Pol1","Pol1 parameters not set");

    return fVal;
}

//______________________________________________________________________________
Double_t TMFunctor::Pol1Exc(Double_t *x, Double_t *par)
{
    if (!fIsFunctorInit)
       Error("TMFunctor::MyFitFunc","Functor not initalized!");

    if (!fReject || (fRejLo == 0. && fRejUp == 0.))
       Error("TMFunctor::Pol1Exc","Rejection from %f to %f not set!",fRejLo,fRejUp);

    if (fReject && x[0] > fRejLo && x[0] < fRejUp) {
      TF1::RejectPoint();
      return 0;
   }
   return par[0] + par[1]*x[0];
}

//______________________________________________________________________________
Double_t TMFunctor::MyFitFuncEmpty(Double_t *x, Double_t *par)
{
    if (!fIsFunctorInit)
       Error("TMFunctor::MyFitFunc","Functor not initalized!");

    Double_t iVal[3];
    Double_t fVal = 0.;

    iVal[0] = hFit[0]->GetBinContent(hFit[0]->FindBin(x[0]));
    iVal[1] = hFit[1]->GetBinContent(hFit[1]->FindBin(x[0]));
    iVal[2] = hFit[2]->GetBinContent(hFit[2]->FindBin(x[0]));

    fVal   = par[0]*iVal[0] + par[1]*iVal[1] + par[2]*iVal[2];

    return fVal;
}

////______________________________________________________________________________
//Double_t TMFunctor::Pol2LowerUpperToMean(Double_t *x, Doubel_t *par)
//{
//    Double_t pol1 = par[0] + par[1]*x[0]+par[2]*x[0]*x[0];
//    Double_t pol2 = par[3] + par[4]*x[0]+par[5]*x[0]*x[0];
//
//    return 0.5*(pol1+pol2);
//}
//
////______________________________________________________________________________
//Double_t TMFunctor::Pol2LowerUpperToMean(Double_t *x, Doubel_t *par)
//{
//    Double_t pol1 = par[0] + par[1]*x[0]+par[2]*x[0]*x[0];
//    Double_t pol2 = par[3] + par[4]*x[0]+par[5]*x[0]*x[0];
//
//    return (pol1-pol2)/2./1.5;
//}

//______________________________________________________________________________
