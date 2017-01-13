/*************************************************************************
 * Author: Manuel Dieterle, 2014
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMCombinedFit                                                                //
//                                                                      //
// Class for determining cuts                                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef myTMCombinedFit
#define myTMCombinedFit

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include "TMFitCollection.h"
#include "TOEnergyThetaCut.h"
#include "TOEnergyThetaData.h"
#include "TNamed.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMCombinedFunctor.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TOLoader.h"

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TEnv.h"
#include "TROOT.h"
#include "TMath.h"


//using namespace std;

class TMCombinedFit : public TNamed
{

protected:

    Int_t nCTBins;            // number of CT bins 

    Bool_t IsW;               // W mode flag

    const Char_t* outName;    // output filename

    TOEnergyThetaData* dBG;   // background factor data
    TOEnergyThetaData* dSig;  // signal factor data
    TOEnergyThetaData* dfBG;  // smoothed background factor data
    TOEnergyThetaData* dfSig; // smoothed signal factor data

    TOEnergyThetaData* dMM;   // missing mass result
    TOEnergyThetaData* dCop;  // coplanarity result
    TOEnergyThetaData* dIM;   // invariant mass result
    TOEnergyThetaData* dfMM;  // smoothed missing mass result
    TOEnergyThetaData* dfCop; // smoothed coplanarity result
    TOEnergyThetaData* dfIM;  // smoothed invariant mass result

    TH1D* hMMbg;              // MM bg
    TH1D* hMMsig;             // MM signal
    TH1D* hCOPbg;             // Cop bg
    TH1D* hCOPsig;            // Cop signal
    TH1D* hIMbg;              // IM bg
    TH1D* hIMsig;             // IM signal
   
    Bool_t kFitCut;           // fit in cut range 
    Bool_t kFixSig;           // fix signal parameter
    Bool_t kFixBG;            // fix bg parameter
    
    Double_t fFitMMLo;        // MM fit minimum range
    Double_t fFitMMUp;        // MM fit maximum range
    Double_t fFitCopLo;       // Cop fit minimum range
    Double_t fFitCopUp;       // Cop fit maximum range
    Double_t fFitIMLo;        // IM fit minimum range
    Double_t fFitIMUp;        // IM fit maximum range

    Double_t fSiglo;          // signal par minimum
    Double_t fSigup;          // signal par maximum
    Double_t fBGlo;           // bg par minimum   
    Double_t fBGup;           // bg par maximum

    const Char_t* fdBG;       // BG factor fit func
    const Char_t* fdSig;      // Sig factor fit func

    Double_t fdSiglo;         // dilfac sig fitting min
    Double_t fdSigup;         // dilfac sig fitting max
    Double_t fdBGlo;          // dilfac fitting min
    Double_t fdBGup;          // dilfac fitting max

    Int_t rebMM;              // MM rebin
    Int_t rebCop;             // Cop rebin
    Int_t rebIM;              // IM rebin

    const Char_t* fRefCS;     // reference cross section

    const Char_t* fTot;       // total data filename
    const Char_t* fBG;        // total bg filename
    const Char_t* fSig;       // total sig filename

    const Char_t* fEffTot;    // total eff filename
    const Char_t* fEffBG;     // bg eff filename
    const Char_t* fEffSig;    // sig eff filename

    const Char_t* fFlux;      // total flux name
    const Char_t* fEff;       // total flux name

    Double_t fTotFluxScale;   // total flux scale factor
    Double_t fBGFluxScale;    // background flux scale factor
    Double_t fSigFluxScale;   // signal flux scale factor

    Double_t fTotNorm;        // normalization factor
    Double_t fBGNorm;         // background normalization factor
    Double_t fSigNorm;        // signal normalization factor

    const Char_t* MMname;     // MM histo name
    const Char_t* Copname;    // Cop histo name
    const Char_t* IMname;     // IM histo name

    const Char_t* fMMcut;     // MM cut filename
    const Char_t* fCopcut;    // Cop cut filename
    const Char_t* fIMcut;     // IM cut filename

    const Char_t* fMM;        // MM cut name
    const Char_t* fCop;       // Cop cut name
    const Char_t* fIM;        // IM cut name

    Bool_t kMMcutW;           // MM cut W flag
    Bool_t kCopcutW;          // Cop cut W flag
    Bool_t kIMcutW;           // IM cut W flag

    Double_t sigMM;           // MM sigma in cut
    Double_t sigCop;          // Cop sigma in cut
    Double_t sigIM;           // IM sigma in cut
    Double_t sigMMuse;        // MM sigma applied
    Double_t sigCopuse;       // Cop sigma applied
    Double_t sigIMuse;        // IM sigma applied

    Double_t fStep;           // step size for fitter
public:

    TMCombinedFit() : TNamed(), nCTBins(0), IsW(0), outName(0), dBG(0), 
                      dSig(0), dfBG(0), dfSig(0), dMM(0), dCop(0), 
                      dIM(0), dfMM(0), dfCop(0), dfIM(0),
                      hMMbg(0), hMMsig(0), hCOPbg(0), hCOPsig(0),
                      hIMbg(0), hIMsig(0), kFitCut(0), kFixSig(0), kFixBG(0),
                      fFitMMLo(0), fFitMMUp(0), fFitCopLo(0), fFitCopUp(0),
                      fFitIMLo(0), fFitIMUp(0), fSiglo(0), fSigup(0),
                      fBGlo(0), fBGup(0), fdBG(0), fdSig(0), 
                      fdSiglo(0), fdSigup(0), fdBGlo(0), fdBGup(0),
                      rebMM(0), rebCop(0), rebIM(0), fRefCS(0),
                      fTot(0), fBG(0), fSig(0), fEffTot(0), fEffBG(0), fEffSig(0),
                      fFlux(0), fEff(0), fTotFluxScale(0), fBGFluxScale(0),
                      fSigFluxScale(0), fTotNorm(0), fBGNorm(0), fSigNorm(0),
                      MMname(0), Copname(0), IMname(0), fMMcut(0), fCopcut(0), 
                      fIMcut(0), fMM(0), fCop(0), fIM(0), kMMcutW(0),
                      kCopcutW(0), kIMcutW(0), sigMM(0), sigCop(0), sigIM(0),
                      sigMMuse(0), sigCopuse(0), sigIMuse(0), fStep(0) {}

    TMCombinedFit(const Char_t* szName, const Char_t* szTitle);
    virtual ~TMCombinedFit();

    void SetIsW(Bool_t k = kTRUE){ IsW = k; }

    void SetMMcut( const Char_t* f, const Char_t* c, Bool_t k = kFALSE ){ fMMcut = f; fMM = c; kMMcutW = k; }
    void SetCopcut( const Char_t* f, const Char_t* c, Bool_t k = kFALSE ){ fCopcut = f; fCop = c; kCopcutW = k; }
    void SetIMcut( const Char_t* f, const Char_t* c, Bool_t k = kFALSE ){ fIMcut = f; fIM = c; kIMcutW = k; }

    void SetMMsigma( Double_t is, Double_t use ){ sigMM = is; sigMMuse = use; } 
    void SetCopsigma( Double_t is, Double_t use ){ sigCop = is; sigCopuse = use; } 
    void SetIMsigma( Double_t is, Double_t use ){ sigIM = is; sigIMuse = use; } 

    void SetFitCutRange( Bool_t k ){ kFitCut = k; }
    void FixSignal( Bool_t k ){ kFixSig = k; }
    void FixBackground( Bool_t k ){ kFixBG = k; }
    void SetMMFitRange( Double_t min, Double_t max ){ fFitMMLo = min; fFitMMUp = max; }
    void SetCopFitRange( Double_t min, Double_t max ){ fFitCopLo = min; fFitCopUp = max; }
    void SetIMFitRange( Double_t min, Double_t max ){ fFitIMLo = min; fFitIMUp = max; }
    void SetStepSize( Double_t f ){ fStep = f; }

    void SetSignalParRange( Double_t min, Double_t max ){ fSiglo = min; fSigup = max; }
    void SetBackgroundParRange( Double_t min, Double_t max ){ fBGlo = min; fBGup = max; }
    void SetDilSignalParRange( Double_t min, Double_t max ){ fdSiglo = min; fdSigup = max; }
    void SetDilBackgroundParRange( Double_t min, Double_t max ){ fdBGlo = min; fdBGup = max; }
    void SetDilFitFuncs( const Char_t* bg, const Char_t* sig ){ fdBG = bg; fdSig = sig; }

    void RebinMMHisto( Int_t r ){ rebMM = r; }
    void RebinCopHisto( Int_t r){ rebCop = r; }
    void RebinIMHisto( Int_t r ){ rebIM = r; }

    void SetReferenceCrossSection( const Char_t* f ){ fRefCS = f; }

    void SetTotalDataFile( const Char_t* f ){ fTot = f; }
    void SetBackgroundDataFile( const Char_t* f ){ fBG = f; }
    void SetSignalDataFile( const Char_t* f ){ fSig = f; }
    void SetTotalEffFile( const Char_t* f ){ fEffTot = f; }
    void SetBackgroundEffFile( const Char_t* f ){ fEffBG = f; }
    void SetSignalEffFile( const Char_t* f ){ fEffSig = f; }

    void SetFluxName( const Char_t* f ){ fFlux = f; }
    void SetEffName( const Char_t* f ){ fEff = f; }

    void ScaleTotalFlux( Double_t f ){ fTotFluxScale = f; }
    void ScaleBackgroundFlux( Double_t f ){ fBGFluxScale = f; }
    void ScaleSignalFlux( Double_t f ){ fSigFluxScale = f; }

    void ScaleTotalHisto( Double_t f ){ fTotNorm = f; }
    void ScaleBackgroundHisto( Double_t f ){ fBGNorm = f; }
    void ScaleSignalHisto( Double_t f ){ fSigNorm = f; }

    void SetMMname( const Char_t* f ){ MMname = f; }
    void SetCopname( const Char_t* f ){ Copname = f; }
    void SetIMname( const Char_t* f ){ IMname = f; }

    void SetNumCosThetaBin( Int_t n ){ nCTBins = n; }

    void DivideFlux(TH2* h, TH1* hFlux);
    void DivideEfficiency(TH2** h, TOEnergyThetaData* eff);
    void Fit();

    Double_t* GetMinMaxRange(TH1* hA, TH1* hB);

    struct GlobalChi2 {
       GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
                    ROOT::Math::IMultiGenFunction & f2,
                    ROOT::Math::IMultiGenFunction & f3) :
          fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3) {}
    
       double operator() (const double *par) const {
    
          return (*fChi2_1)(par) + (*fChi2_2)(par) + (*fChi2_3)(par);
       }
    
       const  ROOT::Math::IMultiGenFunction * fChi2_1;
       const  ROOT::Math::IMultiGenFunction * fChi2_2;
       const  ROOT::Math::IMultiGenFunction * fChi2_3;
    };

    ClassDef(TMCombinedFit, 1)
};

#endif

