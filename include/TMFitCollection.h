/*************************************************************************
 * Author: Manuel Dieterle, 2012
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMFitCollection                                                   //
//                                                                      //
// Class for collecting any kind of TObjects                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef myTMFitCollection
#define myTMFitCollection

#include <cstdlib>
#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TOEnergyThetaCut.h"
#include "TOEnergyThetaData.h"

class TMFitCollection : public TNamed
{

protected:
    Int_t fNBG;                             // size of BG array
    Int_t fNCTbins;                         // # CTbins
    Int_t fNEbins;                          // # Ebins
    Int_t fRep;                             // number of fit repetitions
//    Int_t fIt;                              // applied iteration

    TH2F**   oDatIn;                        // [fNCTbins] input data object
    TH2F**   oSigIn;                        // [fNCTbins] input sig object
    TH2F***  oBGIn;                         // [fNCTbins][fNBG] input sig object
    TH1D***  oDat;                          // [fNCTbins][fNEbins] data object
    TH1D***  oSig;                          // [fNCTbins][fNEbins] sig object
    TH1D**** oBG;                           // [fNCTbins][fNEbins][fNBG] bg object
    TH1D***  oBGTot;                        // [fNCTbins][fNEbins] total bg object
    TH1D***  oTot;                          // [fNCTbins][fNEbins] total object
    TF1***   oFit;                          // [fNCTbins][fNEbins] fit func

    Double_t*    fEnergy;                   // [fNEbins] energy values
    Double_t*    fCTlo;                     // [fNCTbins] lower CT values
    Double_t*    fCTup;                     // [fNCTbins] upper CT values
    Double_t***  fPar;                      // [fNCTbins][fNEbins][fNBG+1] fit parameters
    Double_t***  fParErr;                   // [fNCTbins][fNEbins][fNBG+1] fit parameter Errors
    Double_t** fChiSquare;                  // [fNCTbins][fNEbins] fit chi2
    Double_t fChiMax;                       // maximum chi2 value

    TF1***     fFitCT;                      // [fNEbins][fNBG+1] parameter fit
    TF1***     fFitE;                       // [fNCTbins][fNBG+1] parameter fit error
    TF1** fCutMeanE;                        // [fNCTbins] mean fit func for cut
    TF1** fCutSigmaE;                       // [fNCTbins] sigma fit func for cut
    TF1** fCutLoE;                          // [fNCTbins] lower fit func for cut
    TF1** fCutUpE;                          // [fNCTbins] upper fit func for cut
    TF1** fCutMeanCT;                       // [fNEbins] mean fit func for cut
    TF1** fCutSigmaCT;                      // [fNEbins] sigma fit func for cut
    TF1** fCutLoCT;                         // [fNEbins] lower fit func for cut
    TF1** fCutUpCT;                         // [fNEbins] upper fit func for cut

    TGraphErrors*** gParCT;                  // [fNEbins][fNBG+1] fit parameter graphs for CT
    TGraphErrors*** gParE;                   // [fNCTbins][fNBG+1] fit parameter graphs for E
    TGraphErrors*** gSBCT;                   // [fNEbins][fNBG+1] sb graphs for CT
    TGraphErrors*** gSBE;                    // [fNCTbins][fNBG+1] sb graphs for E
    TGraphErrors**  gCutMeanE;               // [fNCTbins] mean fit graphs for cut
    TGraphErrors**  gCutSigmaE;              // [fNCTbins] sigma fit graphs for cut
    TGraphErrors**  gCutLoE;                 // [fNCTbins] lower fit graphs for cut
    TGraphErrors**  gCutUpE;                 // [fNCTbins] upper fit graphs for cut
    TGraphErrors**  gCutMeanCT;              // [fNEbins] mean fit graphs for cut
    TGraphErrors**  gCutSigmaCT;             // [fNEbins] sigma fit graphs for cut
    TGraphErrors**  gCutLoCT;                // [fNEbins] lower fit graphs for cut
    TGraphErrors**  gCutUpCT;                // [fNEbins] upper fit graphs for cut
    TGraphErrors** gChiSquareE;              // [fNCTbins]
    TGraphErrors** gChiSquareCT;             // [fNEbins]

    Char_t** szPolMeanE;                     // [fNCTbins] fit pol mean cut region
    Char_t** szPolSigmaE;                    // [fNCTbins] fit pol sigma cut region
    Char_t** szPolLoE;                       // [fNCTbins] fit pol lower cut region
    Char_t** szPolUpE;                       // [fNCTbins] fit pol upper cut region
    Char_t* szPolMeanCT;                     // fit pol mean cut region
    Char_t* szPolSigmaCT;                    // fit pol sigma cut region
    Char_t* szPolLoCT;                       // fit pol lower cut region
    Char_t* szPolUpCT;                       // fit pol upper cut region

    Double_t fGmin;                          // gauss fit min range
    Double_t fGmax;                          // gauss fit max range
    Double_t fGMmin;                         // gauss fit mean min 
    Double_t fGMmax;                         // gauss fit mean max
    Double_t fGMtol;                         // gauss fit mean tol
    Double_t fGSmin;                         // gauss fit sigma min
    Double_t fGSmax;                         // gauss fit sigma max
    Double_t fGHmin;                         // gauss fit height min
    Double_t fGHmax;                         // gauss fit height max

    Double_t fMin;                           // fit start value
    Double_t fMax;                           // fit end value
    Double_t* fLo;                           // [fNCTbins] lower fit range
    Double_t* fUp;                           // [fNCTbins] upper fit range
    Double_t fSigma;                         // cut sigma    

    TOEnergyThetaCut* eCut;                  // cut object
    TOEnergyThetaData* eff;                  // efficiency data
    TOEnergyThetaData* sb;                   // signal background ratio
    Double_t fBR;                            // branching ratio
    Double_t fTD;                            // target density
    Double_t fIM;                            // invariant mass dummy value
    const Char_t* sExFunc;                   // ex func histo name
    Bool_t kAll;                             // ex func integration boolean
    TH2F** hIM;                              // [fNCTbins] ExFunc histograms

public:
    TMFitCollection() : TNamed(), fNBG(0), fNCTbins(0), fNEbins(0), fRep(0),// fIt(0),
                        oDatIn(0), oSigIn(0), oBGIn(0), oDat(0), oSig(0),
                        oBG(0), oBGTot(0), oTot(0), oFit(0),
                        fEnergy(0), fCTlo(0), fCTup(0), fPar(0), fParErr(0),
                        fChiSquare(0), fChiMax(0), fFitCT(0), fFitE(0), 
                        fCutMeanE(0), fCutSigmaE(0), fCutLoE(0), fCutUpE(0),
                        fCutMeanCT(0), fCutSigmaCT(0), fCutLoCT(0), fCutUpCT(0),
                        gParCT(0), gParE(0), gSBCT(0), gSBE(0),
                        gCutMeanE(0), gCutSigmaE(0), gCutLoE(0), gCutUpE(0),
                        gCutMeanCT(0), gCutSigmaCT(0), gCutLoCT(0), gCutUpCT(0),
                        gChiSquareE(0), gChiSquareCT(0),
                        szPolMeanE(0), szPolSigmaE(0), szPolLoE(0), szPolUpE(0),
                        szPolMeanCT(0), szPolSigmaCT(0), szPolLoCT(0), szPolUpCT(0),
                        fGmin(0), fGmax(0), fGMmin(0), fGMmax(0), fGMtol(0),
                        fGSmin(0), fGSmax(0), fGHmin(0), fGHmax(0),
                        eCut(0), eff(0), sb(0), fBR(0), fTD(0), fIM(0), sExFunc(0), kAll(0),
                        hIM(0) { }

    TMFitCollection(const Char_t* szName, const Char_t* szTitle);
    virtual ~TMFitCollection();

    void InitObjects();

    void SetDimensions(Int_t c, Int_t e, Int_t b, Int_t r, Int_t f){ fNBG = b; fNCTbins = c; fNEbins  = e; fRep = r; }
    Int_t GetNumCosThetaBins(){ return fNCTbins; }
    Int_t GetNumEnergyBins(){ return fNEbins; }
    Int_t GetNumBackground(){ return fNBG; }
    Int_t GetNumIterations(){ return fRep; }

    void SetGaussianOpt(Double_t f,
                        Double_t fmin, Double_t fmax,
                        Double_t fmmin, Double_t fmmax, Double_t fmtol,
                        Double_t fsmin, Double_t fsmax,
                        Double_t fhmin, Double_t fhmax){ fSigma = f;
                                                         fGmin = fmin; 
                                                         fGmax = fmax;
                                                         fGMmin = fmmin; 
                                                         fGMmax = fmmax;
                                                         fGMtol = fmtol;
                                                         fGSmin = fsmin; 
                                                         fGSmax = fsmax;
                                                         fGHmin = fhmin; 
                                                         fGHmax = fhmax; }

    void AddInputDataObject(Int_t c, TH2F* h){ oDatIn[c] = h; }
    void AddInputSignalObject(Int_t c, TH2F* h){ oSigIn[c] = h; }
    void AddInputBackgroundObject(Int_t c, Int_t b, TH2F* h){ oBGIn[c][b] = h; }
    void AddDataObject(Int_t c, Int_t e, TH1D* h){ oDat[c][e] = h; }
    void AddSignalObject(Int_t c, Int_t e, TH1D* h){ oSig[c][e] = h; }
    void AddBackgroundObject(Int_t c, Int_t e, Int_t b, TH1D* h){ oBG[c][e][b] = h; }
    void AddTotalBackgroundObject(Int_t c, Int_t e, TH1D* h){ oBGTot[c][e] = h; }
    void AddFitObject(Int_t c, Int_t e, TF1* f){ oFit[c][e] = f; }

    void SetTotalHistos(Int_t ct, Int_t e);
    void SetPolynomialParameters(TF1* f, Char_t* s);
    Double_t* GetMaximumAndPositionInRange(TH1D* h);

    void SetEnergy(Int_t e, Double_t f){ fEnergy[e] = f; }
    void SetLowerCosTheta(Int_t c, Double_t f){ fCTlo[c] = f; }
    void SetUpperCosTheta(Int_t c, Double_t f){ fCTup[c] = f; }
    void SetFitParameter(Int_t c, Int_t e, Int_t b, Double_t p){ fPar[c][e][b] = p; }
    void SetFitParError(Int_t c, Int_t e, Int_t b, Double_t p){ fParErr[c][e][b] = p; }
    void SetChiSquare(Int_t c, Int_t e, Double_t f){ fChiSquare[c][e] = f; }
    void SetMaxChiSquare(Double_t f){ fChiMax = f; }

    void SetCTParameterFit(Int_t e, Int_t b, TF1* f){ fFitCT[e][b] = f; }
    void SetEParameterFit(Int_t c, Int_t b, TF1* f){ fFitE[c][b] = f; }

    void SetCTParameterGraph(Int_t e, Int_t b, TGraphErrors* g){ gParCT[e][b] = g; }
    void SetEParameterGraph(Int_t c, Int_t b, TGraphErrors* g){ gParE[c][b] = g; }
    void SetSBCTGraph(Int_t e, Int_t b, TGraphErrors* g){ gSBCT[e][b] = g; }
    void SetSBEGraph(Int_t c, Int_t b, TGraphErrors* g){ gSBE[c][b] = g; }

    void SetMeanCutGraphE(Int_t c, TGraphErrors* g){ gCutMeanE[c] = g; }
    void SetSigmaCutGraphE(Int_t c, TGraphErrors* g){ gCutSigmaE[c] = g; }
    void SetLowerCutGraphE(Int_t c, TGraphErrors* g){ gCutLoE[c] = g; }
    void SetUpperCutGraphE(Int_t c, TGraphErrors* g){ gCutUpE[c] = g; }
    void SetMeanCutGraphCT(Int_t e, TGraphErrors* g){ gCutMeanCT[e] = g; }
    void SetSigmaCutGraphCT(Int_t e, TGraphErrors* g){ gCutSigmaCT[e] = g; }
    void SetLowerCutGraphCT(Int_t e, TGraphErrors* g){ gCutLoCT[e] = g; }
    void SetUpperCutGraphCT(Int_t e, TGraphErrors* g){ gCutUpCT[e] = g; }
    void SetChiSquareGraphE(Int_t c, TGraphErrors* g){ gChiSquareE[c] = g; }
    void SetChiSquareGraphCT(Int_t e, TGraphErrors* g){ gChiSquareCT[e] = g; }

    void SetMeanCutFuncE(Int_t c, TF1* f){ fCutMeanE[c] = f; }
    void SetSigmaCutFuncE(Int_t c, TF1* f){ fCutSigmaE[c] = f; }
    void SetLowerCutFuncE(Int_t c, TF1* f){ fCutLoE[c] = f; }
    void SetUpperCutFuncE(Int_t c, TF1* f){ fCutUpE[c] = f; }
    void SetMeanCutFuncCT(Int_t e, TF1* f){ fCutMeanCT[e] = f; }
    void SetSigmaCutFuncCT(Int_t e, TF1* f){ fCutSigmaCT[e] = f; }
    void SetLowerCutFuncCT(Int_t e, TF1* f){ fCutLoCT[e] = f; }
    void SetUpperCutFuncCT(Int_t e, TF1* f){ fCutUpCT[e] = f; }

    void SetFitLowerBoundary(Double_t v, Char_t* lo);
    void SetFitUpperBoundary(Double_t v, Char_t* up);
    void SetFitPolMeanE(Char_t* s);
    void SetFitPolSigmaE(Char_t* s);
    void SetFitPolLoE(Char_t* s);
    void SetFitPolUpE(Char_t* s);
    void SetFitPolMeanCT(Char_t* s){ szPolMeanCT = s; }
    void SetFitPolSigmaCT(Char_t* s){ szPolSigmaCT = s; }
    void SetFitPolLoCT(Char_t* s){ szPolLoCT = s; }
    void SetFitPolUpCT(Char_t* s){ szPolUpCT = s; }

    void SetExFuncCalcOpt(TOEnergyThetaData* e, 
                          Double_t b, Double_t d, Double_t m,
                          const Char_t* s, Bool_t k){ eff = e;
                                                      fBR = b;
                                                      fTD = d;
                                                      fIM = m; 
                                                      sExFunc = s;
                                                      kAll = k; }

    TH2F* GetInputDataObject(Int_t c){ return oDatIn[c]; }
    TH2F* GetInputSignalObject(Int_t c){ return oSigIn[c]; }
    TH2F* GetInputBackgroundObject(Int_t c, Int_t b){ return oBGIn[c][b]; }
    TH1D* GetDataObject(Int_t c, Int_t e){ return oDat[c][e]; }
    TH1D* GetSignalObject(Int_t c, Int_t e){ return oSig[c][e]; }
    TH1D* GetBackgroundObject(Int_t c, Int_t e, Int_t b){ return oBG[c][e][b]; }
    TH1D* GetTotalBackgroundObject(Int_t c, Int_t e){ return oBGTot[c][e]; }
    TH1D* GetTotalObject(Int_t c, Int_t e){ return oTot[c][e]; }
    TF1* GetFitObject(Int_t c, Int_t e){ return oFit[c][e]; }

    TOEnergyThetaData* GetSBObject(){ return sb; }

    Double_t GetEnergy(Int_t e){ return fEnergy[e]; }
    Double_t GetLowerCosTheta(Int_t c){ return fCTlo[c]; }
    Double_t GetUpperCosTheta(Int_t c){ return fCTup[c]; }
    Int_t GetEnergyBin(Double_t f);
    Int_t GetCosThetaBin(Double_t f);
    Double_t GetFitParameter(Int_t c, Int_t e, Int_t b){ return fPar[c][e][b]; }
    Double_t GetFitParError(Int_t c, Int_t e, Int_t b){ return fParErr[c][e][b]; }
    Double_t GetChiSquare(Int_t c, Int_t e){ return fChiSquare[c][e]; }

    TF1* GetCTParameterFit(Int_t e, Int_t b){ return fFitCT[e][b]; }
    TF1* GetEParameterFit(Int_t c, Int_t b){ return fFitE[c][b]; }
    TF1* GetMeanCutFuncE(Int_t c) { return fCutMeanE[c]; }
    TF1* GetSigmaCutFuncE(Int_t c) { return fCutSigmaE[c]; }
    TF1* GetLowerCutFuncE(Int_t c) { return fCutLoE[c]; }
    TF1* GetUpperCutFuncE(Int_t c) { return fCutUpE[c]; }
    TF1* GetMeanCutFuncCT(Int_t e) { return fCutMeanCT[e]; }
    TF1* GetSigmaCutFuncCT(Int_t e) { return (TF1*)fCutSigmaCT[e]; }
    TF1* GetLowerCutFuncCT(Int_t e) { return (TF1*)fCutLoCT[e]; }
    TF1* GetUpperCutFuncCT(Int_t e) { return (TF1*)fCutUpCT[e]; }

    TGraphErrors* GetCTParameterGraph(Int_t e, Int_t b){ return gParCT[e][b]; }
    TGraphErrors* GetEParameterGraph(Int_t c, Int_t b){ return gParE[c][b]; }
    TGraphErrors* GetSBCTGraph(Int_t e, Int_t b){ return gSBCT[e][b]; }
    TGraphErrors* GetSBEGraph(Int_t c, Int_t b){ return gSBE[c][b]; }
    TGraphErrors* GetMeanCutGraphE(Int_t c) { return gCutMeanE[c]; }
    TGraphErrors* GetSigmaCutGraphE(Int_t c) { return gCutSigmaE[c]; }
    TGraphErrors* GetLowerCutGraphE(Int_t c) { return gCutLoE[c]; }
    TGraphErrors* GetUpperCutGraphE(Int_t c) { return gCutUpE[c]; }
    TGraphErrors* GetMeanCutGraphCT(Int_t e) { return gCutMeanCT[e]; }
    TGraphErrors* GetSigmaCutGraphCT(Int_t e) { return gCutSigmaCT[e]; }
    TGraphErrors* GetLowerCutGraphCT(Int_t e) { return gCutLoCT[e]; }
    TGraphErrors* GetUpperCutGraphCT(Int_t e) { return gCutUpCT[e]; }
    TGraphErrors* GetChiSquareGraphE(Int_t c) { return gChiSquareE[c]; }
    TGraphErrors* GetChiSquareGraphCT(Int_t e) { return gChiSquareCT[e]; }

    void GetCutsFromGaussian();
    TGraphErrors* CleanGraphE(Int_t c, TGraphErrors* g);
    void CleanGraphCT(Int_t e, TGraphErrors* g);
    void CalculateExcitationFunctions(Bool_t kEff);
    TH2** RecalculateExcitationFunctions(TOEnergyThetaData* e, const Char_t* s, Double_t fM, Double_t fB, Double_t fT, Bool_t kA);
    TH2** RecalculateExcFuncErrors(TOEnergyThetaData* e, const Char_t* s, Double_t fM, Double_t fB, Double_t fT, Bool_t kA);

    TOEnergyThetaCut* GetCutObject() { return eCut; }

    Double_t GetFitLowerBoundary(Int_t c){ return fLo[c]; }
    Double_t GetFitUpperBoundary(Int_t c){ return fUp[c]; }
    Char_t* GetFitPolMeanE(Int_t c){ return szPolMeanE[c]; }
    Char_t* GetFitPolSigmaE(Int_t c){ return szPolSigmaE[c]; }
    Char_t* GetFitPolLoE(Int_t c){ return szPolLoE[c]; }
    Char_t* GetFitPolUpE(Int_t c){ return szPolUpE[c]; }

    TH2F* GetExcitationFunction(Int_t ct) { return (TH2F*)hIM[ct]; }

    void Draw(Int_t ct, Int_t e, Double_t rLo = 9999., Double_t rUp = 9999.,
              Bool_t kTit = kTRUE, Bool_t kPrint = kFALSE);
    void DrawChiSquareE(Int_t c);
    void DrawChiSquareCT(Int_t e);
    void DrawInputE(Int_t c, const Char_t* s, const Char_t* t);
    void DrawInputCT(Int_t e, const Char_t* s, const Char_t* t);
    void DrawCuts(Int_t c, const Char_t* s);
    void DrawParametersCT(Int_t e, Int_t bg);
    void DrawParametersE(Int_t ct, Int_t bg);


    ClassDef(TMFitCollection, 1)  // fit object collection
};

#endif
