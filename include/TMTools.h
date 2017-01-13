/*************************************************************************
 * Author: Manuel Dieterle, 2012
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMTools                                                              //
//                                                                      //
// Class containing any necessary tools.                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef myTMTools
#define myTMTools

#include <iostream>
#include <fstream>
#include <stdio.h>
#include "TFile.h"
#include "TKey.h"
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
#include "TStyle.h"
#include "TCutG.h"
#include "TOEnergyThetaCut.h"
#include "TOEnergyThetaData.h"
#include "TMEnergyThetaCut.h"
#include "TMFunctor.h"

namespace TMTools
{
    TH1F* Get1DHisto(TString strFile, TString strOldName, TString strNewName);
    TH2F* Get2DHisto(TString strFile, TString strOldName, TString strNewName);
//    TH3F* Get3DHisto(TString strFile, TString strOldName, TString strNewName);

    void Get2DCollectionSize(TFile* f, const Char_t* szColl, Int_t& ia, Int_t& ib);
    Double_t GetFirstNotEmptyBin(TH1* h);
    Double_t GetLastNotEmptyBin(TH1* h);
    void Draw2DCollection(TFile* f, const Char_t* szColl, Int_t ia, Int_t ib, Int_t iNPadX, Int_t iNPadY, Int_t iNDC, Bool_t kFit = kTRUE, Int_t iRebin = 1);
    void Draw2DCollection(TFile* f, const Char_t* szColl, Int_t ia, Int_t ib, Int_t ic, Int_t iNPadX, Int_t iNPadY, Int_t iNDC, Bool_t kFit = kTRUE, Int_t iRebin = 1);
    void Draw3DCollection(const Char_t* szFile, const Char_t* szColl, Int_t iNPadX, Int_t iNPadY, Int_t iNDC);
    void GetSBfromHistoCollection(const Char_t* szFile, const Char_t* szMM, const Char_t* szColl);
    void GraphToPolygonRight(TFile* f, const Char_t* szGraph, Double_t fLeft, Double_t fBottom, Double_t fRight, Double_t fTop, Double_t fSup = 0.);
    void GraphToPolygonTop(TFile* f, const Char_t* szGraph, Double_t fLeft, Double_t fBottom, Double_t fRight);
    void GetMeanAndSigma(TFile* f, const Char_t* szFunc, Int_t ia, Int_t ib);
    void DrawCutQuality(TFile* f, TFile* fData, const Char_t* szHisto, Int_t iCT, Double_t fSigma, const Char_t* szPol, Double_t fLoLeft, Double_t fLoRight, Double_t fUpLeft = 0., Double_t fUpRight = 0.);
    void DrawCopCutQuality(TFile* f, TFile* fData, const Char_t* szHisto, Int_t iCT, Double_t fSigma, const Char_t* szPol, Double_t fLeft, Double_t fRight);
    void CreateExcitationFunction(TFile* f, const Char_t* szObj, TFile* fEff, Int_t nCTbin, Int_t nTCbin, const Char_t* szHisto, Double_t nBR, Double_t nTD);
    void CorrectEmptyTarget(TFile* fCS, TFile* fCSempty, Double_t scale = 1.);
    TOEnergyThetaData* AverageTOEnergyThetaData(TOEnergyThetaData* e1, TOEnergyThetaData* e2);
    TOEnergyThetaData* AddTOEnergyThetaData(TOEnergyThetaData* e1, TOEnergyThetaData* e2);
    void PrintModelData(TGraphErrors* g, const Char_t* fModel);
    void PrintData(TGraphErrors* g, const Char_t* fModel);
    Bool_t IsObjectType(TObject* o, const Char_t* szType);
    Bool_t CheckObject(const Char_t* s, TFile* f);
//    TGraphErrors** GetLegendreRatios(TOEnergyThetaData* e, TList& list, Int_t pol, 
//                                    Int_t eStart, Int_t eEnd, Double_t left, Double_t right);
}

#endif

