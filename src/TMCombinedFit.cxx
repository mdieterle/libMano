/*************************************************************************
 * Author: Manuel Dieterle, 2012
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMCombinedFit                                                                //
//                                                                      //
// Class for collecting any kind of TObjects                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TMCombinedFit.h"

ClassImp(TMCombinedFit)

//______________________________________________________________________________
TMCombinedFit::TMCombinedFit(const Char_t* szName, const Char_t* szTitle)
    : TNamed(szName, szTitle)
{
    nCTBins = 1;

    outName = szName;

    hMMbg = 0;
    hMMsig = 0;
    hCOPbg = 0;
    hCOPsig = 0;
    hIMbg = 0;
    hIMsig = 0;

    kFitCut = kFALSE;
    kFixSig = kFALSE;
    kFixBG  = kFALSE;
    
    fFitMMLo = 0.;
    fFitMMUp = 0.;
    fFitCopLo = 0.;
    fFitCopUp = 0.;
    fFitIMLo = 0.;
    fFitIMUp = 0.;

    fSiglo = 0.;
    fSigup = 0.;
    fBGlo = 0.;
    fBGup = 0.;

    fdBG  = 0;
    fdSig = 0;

    fdBGlo  = 0.;
    fdBGup  = 0.;
    fdSiglo = 0.;
    fdSigup = 0.;

    rebMM = 0;
    rebCop = 0;
    rebIM = 0;

    fRefCS = "";

    fTot = "";
    fBG = "";
    fSig = "";

    fEffTot = "";
    fEffBG = "";
    fEffSig = "";
    
    fFlux = "";
    fEff = "";
    
    fTotFluxScale = 0.;
    fBGFluxScale = 0.;
    fSigFluxScale = 0.;

    fTotNorm = 0.;
    fBGNorm = 0.;
    fSigNorm = 0.;

    dBG   = 0;
    dSig  = 0;
    dMM   = 0;
    dCop  = 0;
    dIM   = 0;
    dfBG  = 0;
    dfSig = 0;
    dfMM  = 0;
    dfCop = 0;
    dfIM  = 0;

    fStep = 0.01;
}
//______________________________________________________________________________
TMCombinedFit::~TMCombinedFit()
{
}

//_______________________________________________________________________________
void TMCombinedFit::Fit()
{
   Printf("Using %i CT bins",nCTBins);

   // input MM spectra
   TH2** hMM2d = new TH2*[nCTBins];
   TOLoader::LoadObjects(&fTot[0], &MMname[0], &hMM2d, 0, nCTBins-1);
   TH2** hMMbg2d = new TH2*[nCTBins];
   TOLoader::LoadObjects(&fBG[0],  &MMname[0], &hMMbg2d,  0, nCTBins-1);
   TH2** hMMsig2d = new TH2*[nCTBins];
   TOLoader::LoadObjects(&fSig[0], &MMname[0], &hMMsig2d, 0, nCTBins-1);

   // input Cop spectra
   TH2** hCop2d = new TH2*[nCTBins];
   TOLoader::LoadObjects(&fTot[0], &Copname[0], &hCop2d, 0, nCTBins-1);
   TH2** hCopbg2d = new TH2*[nCTBins];
   TOLoader::LoadObjects(&fBG[0],  &Copname[0], &hCopbg2d,  0, nCTBins-1);
   TH2** hCopsig2d = new TH2*[nCTBins];
   TOLoader::LoadObjects(&fSig[0], &Copname[0], &hCopsig2d, 0, nCTBins-1);

   // input IM spectra
   TH2** hIM2d = new TH2*[nCTBins];
   TOLoader::LoadObjects(&fTot[0], &IMname[0], &hIM2d, 0, nCTBins-1);
   TH2** hIMbg2d = new TH2*[nCTBins];
   TOLoader::LoadObjects(&fBG[0],  &IMname[0], &hIMbg2d,  0, nCTBins-1);
   TH2** hIMsig2d = new TH2*[nCTBins];
   TOLoader::LoadObjects(&fSig[0], &IMname[0], &hIMsig2d, 0, nCTBins-1);

   // input flux
   TH1* hFluxTot;
   TOLoader::LoadObject(&fTot[0], &fFlux[0], &hFluxTot);
   TH1* hFluxBG;
   TOLoader::LoadObject(&fBG[0],  &fFlux[0], &hFluxBG);
   TH1* hFluxSig;
   TOLoader::LoadObject(&fSig[0], &fFlux[0], &hFluxSig);

   if (fTotFluxScale)
   {
      hFluxTot->Scale(fTotFluxScale);
      Printf("Scaled dB flux by %f",fTotFluxScale);
   }
   if (fBGFluxScale)
   {
      hFluxBG ->Scale(fBGFluxScale);
      Printf("Scaled C flux by %f",fBGFluxScale);
   }
   if (fSigFluxScale)
   {
      hFluxSig->Scale(fSigFluxScale);
      Printf("Scaled LD2 flux by %f",fSigFluxScale);
   }

   // input eff
   TOEnergyThetaData* effTot;
   TOLoader::LoadObject(&fEffTot[0], &fEff[0], &effTot);
   TOEnergyThetaData* effBG;
   TOLoader::LoadObject(&fEffBG[0],  &fEff[0], &effBG);
   TOEnergyThetaData* effSig;
   TOLoader::LoadObject(&fEffSig[0], &fEff[0], &effSig);

   // reference cross section (only for drawing, not faking :))
   TOEnergyThetaData* refCS;
   TOLoader::LoadObject(&fRefCS[0], "cs", &refCS);

   // prepare dilution factor data
   if (!IsW)
   {
      dBG   = (TOEnergyThetaData*)effTot->Clone("dilFacEg_comb");
      dSig  = (TOEnergyThetaData*)effTot->Clone("sigFacEg_comb");
      dMM   = (TOEnergyThetaData*)effTot->Clone("res_mm");
      dCop  = (TOEnergyThetaData*)effTot->Clone("res_cop");
      dIM   = (TOEnergyThetaData*)effTot->Clone("res_im");

      dfBG  = (TOEnergyThetaData*)effTot->Clone("dilFacEg_smooth_comb");
      dfSig = (TOEnergyThetaData*)effTot->Clone("sigFacEg_smooth_comb");
      dfMM   = (TOEnergyThetaData*)effTot->Clone("res_smooth_mm");
      dfCop  = (TOEnergyThetaData*)effTot->Clone("res_smooth_cop");
      dfIM   = (TOEnergyThetaData*)effTot->Clone("res_smooth_im");
   }
   else
   {
      dBG  = (TOEnergyThetaData*)effTot->Clone("dilFacW_comb");
      dSig = (TOEnergyThetaData*)effTot->Clone("sigFacW_comb");
      dMM  = (TOEnergyThetaData*)effTot->Clone("resw_mm");
      dCop = (TOEnergyThetaData*)effTot->Clone("resw_cop");
      dIM  = (TOEnergyThetaData*)effTot->Clone("resw_im");

      dfBG  = (TOEnergyThetaData*)effTot->Clone("dilFacW_smooth_comb");
      dfSig = (TOEnergyThetaData*)effTot->Clone("sigFacW_smooth_comb");
      dfMM  = (TOEnergyThetaData*)effTot->Clone("resw_smooth_mm");
      dfCop = (TOEnergyThetaData*)effTot->Clone("resw_smooth_cop");
      dfIM  = (TOEnergyThetaData*)effTot->Clone("resw_smooth_im");
   }

   for (Int_t i = 0; i < nCTBins; i++)
   {
      DivideFlux(hMM2d[i],    hFluxTot);
      DivideFlux(hMMbg2d[i],  hFluxBG);
      DivideFlux(hMMsig2d[i], hFluxSig);

      DivideFlux(hCop2d[i],    hFluxTot);
      DivideFlux(hCopbg2d[i],  hFluxBG);
      DivideFlux(hCopsig2d[i], hFluxSig);

      DivideFlux(hIM2d[i],    hFluxTot);
      DivideFlux(hIMbg2d[i],  hFluxBG);
      DivideFlux(hIMsig2d[i], hFluxSig);
    
      if (fTotNorm)
      {
Printf("Scaling dB with 1/%f",1/fTotNorm);
         hMM2d[i]->Scale(fTotNorm);
         hCop2d[i]->Scale(fTotNorm);
         hIM2d[i]->Scale(fTotNorm);
      }

      if (fBGNorm)
      {
Printf("Scaling C with 1/%f",1./fBGNorm);

         hMMbg2d[i]->Scale(fBGNorm);
         hCopbg2d[i]->Scale(fBGNorm);
         hIMbg2d[i]->Scale(fBGNorm);
      }

      if (fSigNorm)
      {
Printf("Scaling LD2 with 1/%f",1./fSigNorm);

         hMMsig2d[i]->Scale(fSigNorm);
         hCopsig2d[i]->Scale(fSigNorm);
         hIMsig2d[i]->Scale(fSigNorm);
      }
   }

   // correct histograms with efficiencies
   DivideEfficiency(hMM2d,     effTot);
   DivideEfficiency(hMMbg2d,   effBG);
   DivideEfficiency(hMMsig2d,  effSig);
   DivideEfficiency(hCop2d,    effTot);
   DivideEfficiency(hCopbg2d,  effBG);
   DivideEfficiency(hCopsig2d, effSig);
   DivideEfficiency(hIM2d,     effTot);
   DivideEfficiency(hIMbg2d,   effBG);
   DivideEfficiency(hIMsig2d,  effSig);

   // get cuts
   Bool_t kUseMM = kFALSE;
   TOEnergyThetaCut* cMM;
   if (strcmp(&fMM[0],"null"))
   {
      TOLoader::LoadObject(&fMMcut[0], &fMM[0], &cMM);
      kUseMM = kTRUE;
   }
   Bool_t kUseCop = kFALSE;
   TOEnergyThetaCut* cCop;
   if (strcmp(&fCop[0],"null"))
   {
      TOLoader::LoadObject(&fCopcut[0], &fCop[0], &cCop);
      kUseCop = kTRUE;
   }
   Bool_t kUseIM = kFALSE;
   TOEnergyThetaCut* cIM;
   if (strcmp(&fIM[0],"null"))
   {
      TOLoader::LoadObject(&fIMcut[0], &fIM[0], &cIM);
      kUseIM = kTRUE;
   }

   // open output file
   TFile* fOut = new TFile(Form("%s",outName),"recreate");

   TCanvas* c2 = 0;
   TLegend* legA = 0;
   TLegend* legB = 0;
   TLine* lMMlo = 0;
   TLine* lMMup = 0;
   TLine* lCoplo = 0;
   TLine* lCopup = 0;
   TLine* lIMlo = 0;
   TLine* lIMup = 0;
   TLine* lMM2lo = 0;
   TLine* lMM2up = 0;
   TLine* lCop2lo = 0;
   TLine* lCop2up = 0;
   TLine* lIM2lo = 0;
   TLine* lIM2up = 0;

   Int_t ipar[2] = {0, 1};

   if (dBG->GetNCosThetaBin() != nCTBins || dBG->GetNEnergyBin() != hMM2d[0]->GetNbinsX())
      Error("TMCombinedFit::Fit","Dilution Factor not of same dimension as input data");

   // Perform Fit
   for (Int_t i = 0; i < nCTBins; i++)
   {
      TGraphErrors* gBG  = new TGraphErrors(hMM2d[i]->GetNbinsX());
      gBG->SetName(Form("carbon_scale_%i",i));
      TGraphErrors* gSig = new TGraphErrors(hMM2d[i]->GetNbinsX());
      gSig->SetName(Form("signal_scale_%i",i));

      for (Int_t j = 0; j < hMM2d[i]->GetNbinsX(); j++)
      {
         Double_t energy    = dBG->GetEnergy(j);
         Double_t Wenergy   = TMath::Sqrt(2.*energy*938.272+938.272*938.272);
         Double_t energyMM  = energy;
         Double_t energyCop = energy;
         Double_t energyIM  = energy;
         if (kMMcutW) energyMM   = Wenergy;
         if (kCopcutW) energyCop = Wenergy;
         if (kIMcutW) energyIM   = Wenergy;

         Double_t CT      = dBG->GetCosThetaBinCenter(i);
         Double_t deltaCT = (Double_t)1./nCTBins;
         Double_t CTlo    = CT-deltaCT; 
         Double_t CTup    = CT+deltaCT;

         Double_t cutMMlo  = 0.;
         Double_t cutMMup  = 0.;

         if (kUseMM)
         {
            cutMMlo = cMM->GetLower(energy,CT);
            cutMMup = cMM->GetUpper(energy,CT);

            if (sigMMuse)
            {
               cutMMlo = cMM->GetMean(energy,CT,sigMM)-sigMMuse*cMM->GetSigma(energy,CT,sigMM);
               cutMMup = cMM->GetMean(energy,CT,sigMM)+sigMMuse*cMM->GetSigma(energy,CT,sigMM);
            }
         }

         Double_t cutCoplo = 0.;
         Double_t cutCopup = 0.;

         if (kUseCop)
         {
            cutCoplo = cCop->GetLower(energy,CT);
            cutCopup = cCop->GetUpper(energy,CT);

            if (sigCopuse)
            {
               cutCoplo = cCop->GetMean(energy,CT,sigCop)-sigCopuse*cCop->GetSigma(energy,CT,sigCop);
               cutCopup = cCop->GetMean(energy,CT,sigCop)+sigCopuse*cCop->GetSigma(energy,CT,sigCop);
            }
         }

         Double_t cutIMlo  = 0.;
         Double_t cutIMup  = 0.;

         if (kUseIM)
         {
            cutIMlo = cIM->GetLower(energy,CT);
            cutIMup = cIM->GetUpper(energy,CT);

            if (sigIMuse)
            {
               cutIMlo = cIM->GetMean(energy,CT,sigIM)-sigIMuse*cIM->GetSigma(energy,CT,sigIM);
               cutIMup = cIM->GetMean(energy,CT,sigIM)+sigIMuse*cIM->GetSigma(energy,CT,sigIM);
            }
         }

         TH1D* hMM  = (TH1D*)hMM2d[i] ->ProjectionY(Form("hMM_%i_%i",i,j), j+1,j+1,"e")->Rebin(rebMM);
         TH1D* hCOP = (TH1D*)hCop2d[i]->ProjectionY(Form("hCop_%i_%i",i,j),j+1,j+1,"e")->Rebin(rebCop);
         TH1D* hIM  = (TH1D*)hIM2d[i] ->ProjectionY(Form("hIM_%i_%i",i,j), j+1,j+1,"e")->Rebin(rebIM);

         hMMbg  = (TH1D*)hMMbg2d[i] ->ProjectionY(Form("hMMbg_%i_%i",i,j), j+1,j+1,"e")->Rebin(rebMM);
         hCOPbg = (TH1D*)hCopbg2d[i]->ProjectionY(Form("hCOPbg_%i_%i",i,j),j+1,j+1,"e")->Rebin(rebCop);
         hIMbg  = (TH1D*)hIMbg2d[i] ->ProjectionY(Form("hIMbg_%i_%i",i,j), j+1,j+1,"e")->Rebin(rebIM);

         hMMsig  = (TH1D*)hMMsig2d[i] ->ProjectionY(Form("hMMsig_%i_%i",i,j), j+1,j+1,"e")->Rebin(rebMM);
         hCOPsig = (TH1D*)hCopsig2d[i]->ProjectionY(Form("hCOPsig_%i_%i",i,j),j+1,j+1,"e")->Rebin(rebCop);
         hIMsig  = (TH1D*)hIMsig2d[i] ->ProjectionY(Form("hIMsig_%i_%i",i,j), j+1,j+1,"e")->Rebin(rebIM);

         TMCombinedFunctor* fCMM  = new TMCombinedFunctor(hMMbg, hMMsig, Form("fCMM_%i_%i",i,j), Form("fCMM_%i_%i",i,j));
         TMCombinedFunctor* fCCop = new TMCombinedFunctor(hCOPbg,hCOPsig,Form("fCCop_%i_%i",i,j),Form("fCCop_%i_%i",i,j));
         TMCombinedFunctor* fCIM  = new TMCombinedFunctor(hIMbg, hIMsig, Form("fCIM_%i_%i",i,j), Form("fCIM_%i_%i",i,j));

         TF1* fMM = 0;
         TF1* fCop = 0;
         TF1* fIM  = 0;

         fMM  = new TF1(Form("fMM_%i_%i",i,j),  fCMM,  &TMCombinedFunctor::MyFitFunc,  fFitMMLo,  fFitMMUp,  2);
         fCop = new TF1(Form("fCop_%i_%i",i,j), fCCop, &TMCombinedFunctor::MyFitFunc,  fFitCopLo, fFitCopUp, 2);
         fIM  = new TF1(Form("fIM_%i_%i",i,j),  fCIM,  &TMCombinedFunctor::MyFitFunc,  fFitIMLo,  fFitIMUp,  2);

         ROOT::Math::WrappedMultiTF1 wfMM( *fMM, 1);
         ROOT::Math::WrappedMultiTF1 wfCop(*fCop,1);
         ROOT::Math::WrappedMultiTF1 wfIM( *fIM, 1);

         ROOT::Fit::DataOptions opt;
         ROOT::Fit::DataRange rangeMM;
         ROOT::Fit::DataRange rangeCop;
         ROOT::Fit::DataRange rangeIM;

         if (kFitCut)
         {
            Printf("Fitting only in cut ranges");
            rangeMM.SetRange(cutMMlo, cutMMup);
            rangeCop.SetRange(cutCoplo, cutCopup);
            rangeIM.SetRange(cutIMlo, cutIMup);
         }
         else
         {
            Printf("Fitting in defined ranges");
            rangeMM.SetRange(fFitMMLo,fFitMMUp);
            rangeCop.SetRange(fFitCopLo,fFitCopUp);
            rangeIM.SetRange(fFitIMLo,fFitIMUp);
         }

         ROOT::Fit::BinData dataMM(opt, rangeMM);
         ROOT::Fit::FillData(dataMM, hMM);

         ROOT::Fit::BinData dataCop(opt, rangeCop);
         ROOT::Fit::FillData(dataCop, hCOP);

         ROOT::Fit::BinData dataIM(opt, rangeIM);
         ROOT::Fit::FillData(dataIM, hIM);

         ROOT::Fit::Chi2Function chi2_MM(dataMM, wfMM);
         ROOT::Fit::Chi2Function chi2_Cop(dataCop, wfCop);
         ROOT::Fit::Chi2Function chi2_IM(dataIM, wfIM);

         GlobalChi2 globalChi2(chi2_MM, chi2_Cop, chi2_IM);

         ROOT::Fit::Fitter fitter;

         const Int_t Npar = 2;
         Double_t par0[Npar] = { 1., 1. };

         fitter.Config().SetParamsSettings(2,par0);
         if (fBGlo && fBGup) fitter.Config().ParSettings(0).SetLimits(fBGlo,   fBGup);
         if (fSiglo && fSigup)fitter.Config().ParSettings(1).SetLimits(fSiglo, fSigup);
         if (kFixBG)  fitter.Config().ParSettings(0).Fix();
         if (kFixSig) fitter.Config().ParSettings(1).Fix();
         fitter.Config().ParSettings(0).SetStepSize(fStep);
         fitter.Config().ParSettings(1).SetStepSize(fStep);

         fitter.Config().MinimizerOptions().SetPrintLevel(0);
         fitter.Config().SetMinimizer("Minuit","Migrad");

         fitter.FitFCN(2,globalChi2,0,dataMM.Size()+dataCop.Size()+dataIM.Size(),true);
         ROOT::Fit::FitResult result = fitter.Result();
         result.Print(std::cout);

         fMM ->SetFitResult( result, ipar);
         fCop->SetFitResult( result, ipar);
         fIM ->SetFitResult( result, ipar);

         dBG ->SetValue(j, i, fMM->GetParameter(0));
         dBG->SetError(j,i,fMM->GetParError(0)); 
         dSig->SetValue(j, i, fMM->GetParameter(1)); 
         dSig->SetError(j,i,fMM->GetParError(1));
         gBG->SetPoint(j,dBG->GetEnergy(j),fMM->GetParameter(0));
         gBG->SetPointError(j,0.,fMM->GetParError(0));

         gSig->SetPoint(j,dSig->GetEnergy(j),fMM->GetParameter(1));
         gSig->SetPointError(j,0.,fMM->GetParError(1));

         hMMbg ->Scale(fMM->GetParameter(0));///hMM->GetMaximum());
         hMMsig->Scale(fMM->GetParameter(1));///hMM->GetMaximum());
//         hMM   ->Scale(1./hMM->GetMaximum());
         TH1D* hMMtot = (TH1D*)hMMbg->Clone(Form("hMMtot_%i_%i",i,j));
         hMMtot->Add(hMMsig);

         hCOPbg ->Scale(fCop->GetParameter(0));///hCOP->GetMaximum());
         hCOPsig->Scale(fCop->GetParameter(1));///hCOP->GetMaximum());
//         hCOP   ->Scale(1./hCOP->GetMaximum());
         TH1D* hCOPtot = (TH1D*)hCOPbg->Clone(Form("hCOPtot_%i_%i",i,j));
         hCOPtot->Add(hCOPsig);

         hIMbg ->Scale(fIM->GetParameter(0));///hIM->GetMaximum());
         hIMsig->Scale(fIM->GetParameter(1));///hIM->GetMaximum());
//         hIM   ->Scale(1./hIM->GetMaximum());
         TH1D* hIMtot = (TH1D*)hIMbg->Clone(Form("hIMtot_%i_%i",i,j));
         hIMtot->Add(hIMsig);

         hMMtot->SetLineColor(kRed);
         hMMbg->SetLineColor(kBlue);
         hMMsig->SetLineColor(kGreen);
         hCOPtot->SetLineColor(kRed);
         hCOPbg->SetLineColor(kBlue);
         hCOPsig->SetLineColor(kGreen);
         hIMtot->SetLineColor(kRed);
         hIMbg->SetLineColor(kBlue);
         hIMsig->SetLineColor(kGreen);

         hMMtot->SetLineWidth(2);
         hMMbg->SetLineWidth(2);
         hMMsig->SetLineWidth(2);
         hCOPtot->SetLineWidth(2);
         hCOPbg->SetLineWidth(2);
         hCOPsig->SetLineWidth(2);
         hIMtot->SetLineWidth(2);
         hIMbg->SetLineWidth(2);
         hIMsig->SetLineWidth(2);

         TH1D* hMMDiff = (TH1D*)hMM->Clone(Form("hMMDiff_%i_%i",i,j));
         hMMDiff->Add(hMMbg,-1.);

         TH1D* hCOPDiff = (TH1D*)hCOP->Clone(Form("hCopDiff_%i_%i",i,j));
         hCOPDiff->Add(hCOPbg,-1.);

         TH1D* hIMDiff = (TH1D*)hIM->Clone(Form("hIMDiff_%i_%i",i,j));
         hIMDiff->Add(hIMbg,-1.);

         hMMDiff->SetMarkerColor(15);
         hCOPDiff->SetMarkerColor(15);
         hIMDiff->SetMarkerColor(15);

         hMMDiff->SetLineWidth(2);
         hCOPDiff->SetLineWidth(2);
         hIMDiff->SetLineWidth(2);

         // write histos for plotting
         hMM->Write();
         hMMtot->Write();
         hMMbg->Write();
         hMMsig->Write();
         hMMDiff->Write();
         hIM->Write();
         hIMtot->Write();
         hIMbg->Write();
         hIMsig->Write();
         hIMDiff->Write();
         hCOP->Write();
         hCOPtot->Write();
         hCOPbg->Write();
         hCOPsig->Write();
         hCOPDiff->Write();

         legA = new TLegend(0.1,0.9,0.9,0.1);
         legA->SetTextFont(41);
         legA->SetFillColor(0);
         legA->AddEntry(hMM,    "dButanol","l");
         legA->AddEntry(hMMsig, Form("LD2*%f",fMM->GetParameter(1)),"l");
         legA->AddEntry(hMMbg,  Form("Carbon*%f",fMM->GetParameter(0)),"l");
         legA->AddEntry(hMMtot, "LD2+Carbon","l");

         legB = new TLegend(0.1,0.9,0.9,0.1);
         legB->SetTextFont(41);
         legB->SetFillColor(0);
         legB->AddEntry(hMMsig, "LD2","l");
         legB->AddEntry(hMMDiff, "dButanol-Carbon","l");

         c2 = new TCanvas(Form("fitres_%i_%i",i,j),Form("fitres_%i_%i",i,j),1300, 1000);
         c2->Divide(4, 2, 0.0001, 0.0001);

         Double_t* rMM   = GetMinMaxRange(hMM,  hMMtot);
         Double_t* rCop  = GetMinMaxRange(hCOP, hCOPtot);
         Double_t* rIM   = GetMinMaxRange(hIM,  hIMtot);
         Double_t* r2MM  = GetMinMaxRange(hMMDiff,  hMMsig);
         Double_t* r2Cop = GetMinMaxRange(hCOPDiff, hCOPsig);
         Double_t* r2IM  = GetMinMaxRange(hIMDiff,  hIMsig);

         Double_t chi     = result.Chi2()/result.Ndf();

         Double_t intMM    = 0.;

         if (kUseMM)
         {
            Double_t errMM = 0.;
            intMM = hMMDiff->IntegralAndError(hMMDiff->FindBin(cutMMlo),hMMDiff->FindBin(cutMMup),errMM);
            dMM->SetValue(j,i,intMM);
            dMM->SetError(j,i,errMM);

//            intMM = hMM->Integral(hMM->FindBin(cutMMlo),hMM->FindBin(cutMMup))
//                    - hMMbg->Integral(hMMbg->FindBin(cutMMlo),hMMbg->FindBin(cutMMup));

            lMMlo   = new TLine(cutMMlo, rMM[0],
                                cutMMlo, rMM[1]);
            lMMup   = new TLine(cutMMup, rMM[0],
                                cutMMup, rMM[1]);
            lMM2lo  = new TLine(cutMMlo, r2MM[0],
                                cutMMlo, r2MM[1]);
            lMM2up  = new TLine(cutMMup, r2MM[0],
                                cutMMup, r2MM[1]);

            lMMlo->SetLineWidth(2);
            lMMlo->SetLineStyle(2);
            lMMup->SetLineWidth(2);
            lMMup->SetLineStyle(2);
            lMM2lo->SetLineWidth(2);
            lMM2lo->SetLineStyle(2);
            lMM2up->SetLineWidth(2);
            lMM2up->SetLineStyle(2);
         }

         Double_t intCop = 0.;

         if (kUseCop)
         {
            Double_t errCop = 0.;
            intCop = hCOPDiff->IntegralAndError(hCOPDiff->FindBin(cutCoplo),hCOPDiff->FindBin(cutCopup),errCop);
            dCop->SetValue(j,i,intCop);
            dCop->SetError(j,i,errCop);

//            intCop = hCOP->Integral(hCOP->FindBin(cutCoplo),hCOP->FindBin(cutCopup))
//                    - hCOPbg->Integral(hCOPbg->FindBin(cutCoplo),hCOPbg->FindBin(cutCopup));

            lCoplo   = new TLine(cutCoplo, rCop[0],
                                cutCoplo, rCop[1]);
            lCopup   = new TLine(cutCopup, rCop[0],
                                cutCopup, rCop[1]);
            lCop2lo  = new TLine(cutCoplo, r2Cop[0],
                                cutCoplo, r2Cop[1]);
            lCop2up  = new TLine(cutCopup, r2Cop[0],
                                cutCopup, r2Cop[1]);

            lCoplo->SetLineWidth(2);
            lCoplo->SetLineStyle(2);
            lCopup->SetLineWidth(2);
            lCopup->SetLineStyle(2);
            lCop2lo->SetLineWidth(2);
            lCop2lo->SetLineStyle(2);
            lCop2up->SetLineWidth(2);
            lCop2up->SetLineStyle(2);
         }

         Double_t intIM = 0.;

         if (kUseIM)
         {
            Double_t errIM = 0.;
            intIM = hIMDiff->IntegralAndError(hIMDiff->FindBin(cutIMlo),hIMDiff->FindBin(cutIMup),errIM);
            dIM->SetValue(j,i,intIM);
            dIM->SetError(j,i,errIM);

//            intIM = hIM->Integral(hIM->FindBin(cutIMlo),hIM->FindBin(cutIMup))
//                    - hIMbg->Integral(hIMbg->FindBin(cutIMlo),hIMbg->FindBin(cutIMup));

            lIMlo   = new TLine(cutIMlo, rIM[0],
                                cutIMlo, rIM[1]);
            lIMup   = new TLine(cutIMup, rIM[0],
                                cutIMup, rIM[1]);
            lIM2lo  = new TLine(cutIMlo, r2IM[0],
                                cutIMlo, r2IM[1]);
            lIM2up  = new TLine(cutIMup, r2IM[0],
                                cutIMup, r2IM[1]);

            lIMlo->SetLineWidth(2);
            lIMlo->SetLineStyle(2);
            lIMup->SetLineWidth(2);
            lIMup->SetLineStyle(2);
            lIM2lo->SetLineWidth(2);
            lIM2lo->SetLineStyle(2);
            lIM2up->SetLineWidth(2);
            lIM2up->SetLineStyle(2);
         }

         Char_t tit[512];
         sprintf(tit,"E_{#gamma}=%5.1fMeV | %2.1f < cos(#theta) < %2.1f | Chi2/ndf=%f",energy,CTlo,CTup,chi);

         c2->cd(1);
         hMM->GetXaxis()->SetRangeUser(fFitMMLo,fFitMMUp);
         hMM->GetYaxis()->SetRangeUser(rMM[0],rMM[1]);
         hMM->GetXaxis()->SetTitleSize(1.5);
         hMM->SetTitle(tit);
         hMM->Draw("e1");
         hMMbg->Draw("hist same");
         hMMsig->Draw("hist same");
         hMMtot->Draw("hist same");
         if (kUseMM)
         {
            lMMlo->Draw("c same");
            lMMup->Draw("c same");
         }

         c2->cd(2);
         hCOP->GetXaxis()->SetRangeUser(fFitCopLo,fFitCopUp);
         hCOP->GetYaxis()->SetRangeUser(rCop[0],rCop[1]);
         hCOP->GetXaxis()->SetTitleSize(1.5);
         hCOP->SetTitle(tit);
         hCOP->Draw("e1");
         hCOPbg->Draw("hist same");
         hCOPsig->Draw("hist same");
         hCOPtot->Draw("hist same");
         if (kUseCop)
         {
            lCoplo->Draw("c same");
            lCopup->Draw("c same");
         }

         c2->cd(3);
         hIM->GetXaxis()->SetRangeUser(fFitIMLo,fFitIMUp);
         hIM->GetYaxis()->SetRangeUser(rIM[0],rIM[1]);
         hIM->GetXaxis()->SetTitleSize(1.5);
         hIM->SetTitle(tit);
         hIM->Draw("e1");
         hIMbg->Draw("hist same");
         hIMsig->Draw("hist same");
         hIMtot->Draw("hist same");
         if (kUseIM)
         {
            lIMlo->Draw("c same");
            lIMup->Draw("c same");
         }

         c2->cd(4);
         legA->Draw();

         c2->cd(5);
         if (kUseMM) hMMsig->SetTitle(Form("int=%4.3e",intMM));
         else hMMsig->SetTitle("");
         hMMsig->GetXaxis()->SetTitleSize(1.5);
         hMMsig->GetXaxis()->SetRangeUser(fFitMMLo,fFitMMUp);
         hMMsig->GetYaxis()->SetRangeUser(r2MM[0],r2MM[1]);
         hMMsig->Draw("e1");
         hMMDiff->Draw("hist same");
         if (kUseMM)
         {
            lMM2lo->Draw("c same");
            lMM2up->Draw("c same");
         }

         c2->cd(6);
         if (kUseCop) hCOPsig->SetTitle(Form("int=%4.3e",intCop));
         else hCOPsig->SetTitle("");
         hCOPsig->GetXaxis()->SetTitleSize(1.5);
         hCOPsig->GetXaxis()->SetRangeUser(fFitCopLo,fFitCopUp);
         hCOPsig->GetYaxis()->SetRangeUser(r2Cop[0],r2Cop[1]);
         hCOPsig->Draw("e1");
         hCOPDiff->Draw("hist same");
         if (kUseCop)
         {
            lCop2lo->Draw("c same");
            lCop2up->Draw("c same");
         }

         c2->cd(7);
         if (kUseIM) hIMsig->SetTitle(Form("int=%4.3e",intIM));
         else hIMsig->SetTitle("");
         hIMsig->GetXaxis()->SetTitleSize(1.5);
         hIMsig->GetXaxis()->SetRangeUser(fFitIMLo,fFitIMUp);
         hIMsig->GetYaxis()->SetRangeUser(r2IM[0],r2IM[1]);
         hIMsig->Draw("e1");
         hIMDiff->Draw("hist same");
         if (kUseIM)
         {
            lIM2lo->Draw("c same");
            lIM2up->Draw("c same");
         }

         c2->cd(8);
         legB->Draw();

         fOut->cd();
         c2->Write();
      }

      TF1* fBG  = new TF1(Form("fBG_%i",i), &fdBG[0],dBG->GetEnergy(0), dBG->GetEnergy(dBG->GetNEnergyBin()-1));
      TF1* fSig = new TF1(Form("fSig_%i",i),&fdSig[0],dSig->GetEnergy(0),dSig->GetEnergy(dSig->GetNEnergyBin()-1));
      gBG->Fit(fBG,"+R","",fdBGlo,  fdBGup);
      gSig->Fit(fSig,"+R","",fdSiglo, fdSigup);

      for (Int_t j = 0; j < dfBG->GetNEnergyBin(); j++)
      {
         dfBG->SetValue(j,i,fBG->Eval(dfBG->GetEnergy(j)));
         dfBG->SetError(j,i,dBG->GetError(j,i));
         dfSig->SetValue(j,i,fSig->Eval(dfSig->GetEnergy(j)));
         dfSig->SetError(j,i,dSig->GetError(j,i));
      }

      gBG->Write();
      gSig->Write();

      delete gBG;
      delete gSig;
      delete fBG;
      delete fSig;
   }

   dMM->GetData()->Scale(nCTBins/4./TMath::Pi());
   dCop->GetData()->Scale(nCTBins/4./TMath::Pi());
   dIM->GetData()->Scale(nCTBins/4./TMath::Pi());

   dMM->GetData()->Scale(1.e6);
   dCop->GetData()->Scale(1.e6);
   dIM->GetData()->Scale(1.e6);

   TCanvas* c3 = new TCanvas("c3","c3",1300,1000);
   TGraphErrors* gRef = refCS->CreateTotalIntegralSum();
   TGraphErrors* gMM  = dMM->CreateTotalIntegralSum();
   TGraphErrors* gCop = dCop->CreateTotalIntegralSum();
   TGraphErrors* gIM  = dIM->CreateTotalIntegralSum();
   gRef->SetMarkerStyle(20);
   gRef->SetMarkerColor(kBlack);
   gRef->SetMarkerSize(1.5);
   gMM->SetMarkerStyle(20);
   gMM->SetMarkerColor(kBlue);
   gMM->SetMarkerSize(1.5);
   gCop->SetMarkerStyle(20);
   gCop->SetMarkerColor(kGreen);
   gMM->SetMarkerSize(1.5);
   gIM->SetMarkerStyle(20);
   gIM->SetMarkerColor(kRed);
   gIM->SetMarkerSize(1.5);

   TLegend* csl = new TLegend(0.7,0.7,0.85,0.85);
   csl->SetBorderSize(0);
   csl->SetFillColor(0);
   csl->SetTextSize(0.045);
   csl->SetTextFont(42);
   csl->AddEntry(gRef,"ref CS","p");
   csl->AddEntry(gMM,"MM","p");
   csl->AddEntry(gCop,"Cop","p");
   csl->AddEntry(gIM,"IM","p");

   c3->cd();
   gRef->Draw("ap");
   gRef->GetYaxis()->SetRangeUser(-0.5,1.5*TMath::MaxElement(gRef->GetN(),gRef->GetY()));//1.5*gRef->GetMaximum());
   gMM->Draw("psame");
   gCop->Draw("psame");
   gIM->Draw("psame");
   csl->Draw("same");

   c3->Write();

   for (Int_t i = 0; i < dBG->GetNEnergyBin(); i++)
      for (Int_t j = 0; j < dBG->GetNCosThetaBin(); j++)
         Printf("Carbon Factor(corr):%f(%f) LD2 Factor(corr):%f(%f)",dBG->GetValue(i,j),dfBG->GetValue(i,j),
                                                                     dSig->GetValue(i,j),dfSig->GetValue(i,j));

   delete c2;
   delete c3;
   delete csl;
   delete legA;
   delete legB;

   // repeat and apply fitted factors
   for (Int_t i = 0; i < nCTBins; i++)
   {
      for (Int_t j = 0; j < hMM2d[i]->GetNbinsX(); j++)
      {
         Double_t energy    = dBG->GetEnergy(j);
         Double_t Wenergy   = TMath::Sqrt(2.*energy*938.272+938.272*938.272);
         Double_t energyMM  = energy;
         Double_t energyCop = energy;
         Double_t energyIM  = energy;
         if (kMMcutW) energyMM   = Wenergy;
         if (kCopcutW) energyCop = Wenergy;
         if (kIMcutW) energyIM   = Wenergy;

         Double_t CT      = dBG->GetCosThetaBinCenter(i);
         Double_t deltaCT = (Double_t)1./nCTBins;
         Double_t CTlo    = CT-deltaCT; 
         Double_t CTup    = CT+deltaCT;

         Double_t cutMMlo  = 0.;
         Double_t cutMMup  = 0.;

         if (kUseMM)
         {
            cutMMlo = cMM->GetLower(energy,CT);
            cutMMup = cMM->GetUpper(energy,CT);

            if (sigMMuse)
            {
               cutMMlo = cMM->GetMean(energy,CT,sigMM)-sigMMuse*cMM->GetSigma(energy,CT,sigMM);
               cutMMup = cMM->GetMean(energy,CT,sigMM)+sigMMuse*cMM->GetSigma(energy,CT,sigMM);
            }
         }

         Double_t cutCoplo = 0.;
         Double_t cutCopup = 0.;

         if (kUseCop)
         {
            cutCoplo = cCop->GetLower(energy,CT);
            cutCopup = cCop->GetUpper(energy,CT);

            if (sigCopuse)
            {
               cutCoplo = cCop->GetMean(energy,CT,sigCop)-sigCopuse*cCop->GetSigma(energy,CT,sigCop);
               cutCopup = cCop->GetMean(energy,CT,sigCop)+sigCopuse*cCop->GetSigma(energy,CT,sigCop);
            }
         }

         Double_t cutIMlo  = 0.;
         Double_t cutIMup  = 0.;

         if (kUseIM)
         {
            cutIMlo = cIM->GetLower(energy,CT);
            cutIMup = cIM->GetUpper(energy,CT);

            if (sigIMuse)
            {
               cutIMlo = cIM->GetMean(energy,CT,sigIM)-sigIMuse*cIM->GetSigma(energy,CT,sigIM);
               cutIMup = cIM->GetMean(energy,CT,sigIM)+sigIMuse*cIM->GetSigma(energy,CT,sigIM);
            }
         }

         TH1D* hMM  = (TH1D*)hMM2d[i] ->ProjectionY(Form("hfMM_%i_%i",i,j), j+1,j+1,"e")->Rebin(rebMM);
         TH1D* hCOP = (TH1D*)hCop2d[i]->ProjectionY(Form("hfCop_%i_%i",i,j),j+1,j+1,"e")->Rebin(rebCop);
         TH1D* hIM  = (TH1D*)hIM2d[i] ->ProjectionY(Form("hfIM_%i_%i",i,j), j+1,j+1,"e")->Rebin(rebIM);

         hMMbg  = (TH1D*)hMMbg2d[i] ->ProjectionY(Form("hfMMbg_%i_%i",i,j), j+1,j+1,"e")->Rebin(rebMM);
         hCOPbg = (TH1D*)hCopbg2d[i]->ProjectionY(Form("hfCOPbg_%i_%i",i,j),j+1,j+1,"e")->Rebin(rebCop);
         hIMbg  = (TH1D*)hIMbg2d[i] ->ProjectionY(Form("hfIMbg_%i_%i",i,j), j+1,j+1,"e")->Rebin(rebIM);

         hMMsig  = (TH1D*)hMMsig2d[i] ->ProjectionY(Form("hfMMsig_%i_%i",i,j), j+1,j+1,"e")->Rebin(rebMM);
         hCOPsig = (TH1D*)hCopsig2d[i]->ProjectionY(Form("hfCOPsig_%i_%i",i,j),j+1,j+1,"e")->Rebin(rebCop);
         hIMsig  = (TH1D*)hIMsig2d[i] ->ProjectionY(Form("hfIMsig_%i_%i",i,j), j+1,j+1,"e")->Rebin(rebIM);

         hMMbg ->Scale(dfBG->GetValue(j,i));///hMM->GetMaximum());
         hMMsig->Scale(dfSig->GetValue(j,i));///hMM->GetMaximum());
         TH1D* hMMtot = (TH1D*)hMMbg->Clone(Form("hfMMtot_%i_%i",i,j));
         hMMtot->Add(hMMsig);

         hCOPbg ->Scale(dfBG->GetValue(j,i));///hCOP->GetMaximum());
         hCOPsig->Scale(dfSig->GetValue(j,i));///hCOP->GetMaximum());
         TH1D* hCOPtot = (TH1D*)hCOPbg->Clone(Form("hfCOPtot_%i_%i",i,j));
         hCOPtot->Add(hCOPsig);

         hIMbg ->Scale(dfBG->GetValue(j,i));///hIM->GetMaximum());
         hIMsig->Scale(dfSig->GetValue(j,i));///hIM->GetMaximum());
         TH1D* hIMtot = (TH1D*)hIMbg->Clone(Form("hfIMtot_%i_%i",i,j));
         hIMtot->Add(hIMsig);

         hMMtot->SetLineColor(kRed);
         hMMbg->SetLineColor(kBlue);
         hMMsig->SetLineColor(kGreen);
         hCOPtot->SetLineColor(kRed);
         hCOPbg->SetLineColor(kBlue);
         hCOPsig->SetLineColor(kGreen);
         hIMtot->SetLineColor(kRed);
         hIMbg->SetLineColor(kBlue);
         hIMsig->SetLineColor(kGreen);

         hMMtot->SetLineWidth(2);
         hMMbg->SetLineWidth(2);
         hMMsig->SetLineWidth(2);
         hCOPtot->SetLineWidth(2);
         hCOPbg->SetLineWidth(2);
         hCOPsig->SetLineWidth(2);
         hIMtot->SetLineWidth(2);
         hIMbg->SetLineWidth(2);
         hIMsig->SetLineWidth(2);

         TH1D* hMMDiff = (TH1D*)hMM->Clone(Form("hfMMDiff_%i_%i",i,j));
         hMMDiff->Add(hMMbg,-1.);

         TH1D* hCOPDiff = (TH1D*)hCOP->Clone(Form("hfCopDiff_%i_%i",i,j));
         hCOPDiff->Add(hCOPbg,-1.);

         TH1D* hIMDiff = (TH1D*)hIM->Clone(Form("hfIMDiff_%i_%i",i,j));
         hIMDiff->Add(hIMbg,-1.);

         hMMDiff->SetMarkerColor(15);
         hCOPDiff->SetMarkerColor(15);
         hIMDiff->SetMarkerColor(15);

         hMMDiff->SetLineWidth(2);
         hCOPDiff->SetLineWidth(2);
         hIMDiff->SetLineWidth(2);

         legA = new TLegend(0.1,0.9,0.9,0.1);
         legA->SetTextFont(41);
         legA->SetFillColor(0);
         legA->AddEntry(hMM,    "dButanol","l");
         legA->AddEntry(hMMsig, Form("LD2*%f",dfSig->GetValue(j,i)),"l");
         legA->AddEntry(hMMbg,  Form("Carbon*%f",dfBG->GetValue(j,i)),"l");
         legA->AddEntry(hMMtot, "LD2+Carbon","l");

         legB = new TLegend(0.1,0.9,0.9,0.1);
         legB->SetTextFont(41);
         legB->SetFillColor(0);
         legB->AddEntry(hMMsig, "LD2","l");
         legB->AddEntry(hMMDiff, "dButanol-Carbon","l");

         c2 = new TCanvas(Form("fitres_smooth_%i_%i",i,j),Form("fitres_smooth_%i_%i",i,j),1300, 1000);
         c2->Divide(4, 2, 0.0001, 0.0001);

         Double_t* rMM   = GetMinMaxRange(hMM,  hMMtot);
         Double_t* rCop  = GetMinMaxRange(hCOP, hCOPtot);
         Double_t* rIM   = GetMinMaxRange(hIM,  hIMtot);
         Double_t* r2MM  = GetMinMaxRange(hMMDiff,  hMMsig);
         Double_t* r2Cop = GetMinMaxRange(hCOPDiff, hCOPsig);
         Double_t* r2IM  = GetMinMaxRange(hIMDiff,  hIMsig);

         Double_t intMM    = 0.;

         if (kUseMM)
         {
            Double_t errMM = 0.;
            intMM = hMMDiff->IntegralAndError(hMMDiff->FindBin(cutMMlo),hMMDiff->FindBin(cutMMup),errMM);
            dfMM->SetValue(j,i,intMM);
            dfMM->SetError(j,i,errMM);

//            intMM = hMM->Integral(hMM->FindBin(cutMMlo),hMM->FindBin(cutMMup))
//                    - hMMbg->Integral(hMMbg->FindBin(cutMMlo),hMMbg->FindBin(cutMMup));

            lMMlo   = new TLine(cutMMlo, rMM[0],
                                cutMMlo, rMM[1]);
            lMMup   = new TLine(cutMMup, rMM[0],
                                cutMMup, rMM[1]);
            lMM2lo  = new TLine(cutMMlo, r2MM[0],
                                cutMMlo, r2MM[1]);
            lMM2up  = new TLine(cutMMup, r2MM[0],
                                cutMMup, r2MM[1]);

            lMMlo->SetLineWidth(2);
            lMMlo->SetLineStyle(2);
            lMMup->SetLineWidth(2);
            lMMup->SetLineStyle(2);
            lMM2lo->SetLineWidth(2);
            lMM2lo->SetLineStyle(2);
            lMM2up->SetLineWidth(2);
            lMM2up->SetLineStyle(2);
         }

         Double_t intCop = 0.;

         if (kUseCop)
         {
            Double_t errCop = 0.;
            intCop = hCOPDiff->IntegralAndError(hCOPDiff->FindBin(cutCoplo),hCOPDiff->FindBin(cutCopup),errCop);
            dfCop->SetValue(j,i,intCop);
            dfCop->SetError(j,i,errCop);

//            intCop = hCOP->Integral(hCOP->FindBin(cutCoplo),hCOP->FindBin(cutCopup))
//                    - hCOPbg->Integral(hCOPbg->FindBin(cutCoplo),hCOPbg->FindBin(cutCopup));

            lCoplo   = new TLine(cutCoplo, rCop[0],
                                cutCoplo, rCop[1]);
            lCopup   = new TLine(cutCopup, rCop[0],
                                cutCopup, rCop[1]);
            lCop2lo  = new TLine(cutCoplo, r2Cop[0],
                                cutCoplo, r2Cop[1]);
            lCop2up  = new TLine(cutCopup, r2Cop[0],
                                cutCopup, r2Cop[1]);

            lCoplo->SetLineWidth(2);
            lCoplo->SetLineStyle(2);
            lCopup->SetLineWidth(2);
            lCopup->SetLineStyle(2);
            lCop2lo->SetLineWidth(2);
            lCop2lo->SetLineStyle(2);
            lCop2up->SetLineWidth(2);
            lCop2up->SetLineStyle(2);
         }

         Double_t intIM = 0.;

         if (kUseIM)
         {
            Double_t errIM = 0.;
            intIM = hIMDiff->IntegralAndError(hIMDiff->FindBin(cutIMlo),hIMDiff->FindBin(cutIMup),errIM);
            dfIM->SetValue(j,i,intIM);
            dfIM->SetError(j,i,errIM);

//            intIM = hIM->Integral(hIM->FindBin(cutIMlo),hIM->FindBin(cutIMup))
//                    - hIMbg->Integral(hIMbg->FindBin(cutIMlo),hIMbg->FindBin(cutIMup));

            lIMlo   = new TLine(cutIMlo, rIM[0],
                                cutIMlo, rIM[1]);
            lIMup   = new TLine(cutIMup, rIM[0],
                                cutIMup, rIM[1]);
            lIM2lo  = new TLine(cutIMlo, r2IM[0],
                                cutIMlo, r2IM[1]);
            lIM2up  = new TLine(cutIMup, r2IM[0],
                                cutIMup, r2IM[1]);

            lIMlo->SetLineWidth(2);
            lIMlo->SetLineStyle(2);
            lIMup->SetLineWidth(2);
            lIMup->SetLineStyle(2);
            lIM2lo->SetLineWidth(2);
            lIM2lo->SetLineStyle(2);
            lIM2up->SetLineWidth(2);
            lIM2up->SetLineStyle(2);
         }

         Char_t tit[512];
         sprintf(tit,"E_{#gamma}=%5.1fMeV | %2.1f < cos(#theta) < %2.1f",energy,CTlo,CTup);

         c2->cd(1);
         hMM->GetXaxis()->SetRangeUser(fFitMMLo,fFitMMUp);
         hMM->GetYaxis()->SetRangeUser(rMM[0],rMM[1]);
         hMM->GetXaxis()->SetTitleSize(1.5);
         hMM->SetTitle(tit);
         hMM->Draw("e1");
         hMMbg->Draw("hist same");
         hMMsig->Draw("hist same");
         hMMtot->Draw("hist same");
         if (kUseMM)
         {
            lMMlo->Draw("c same");
            lMMup->Draw("c same");
         }

         c2->cd(2);
         hCOP->GetXaxis()->SetRangeUser(fFitCopLo,fFitCopUp);
         hCOP->GetYaxis()->SetRangeUser(rCop[0],rCop[1]);
         hCOP->GetXaxis()->SetTitleSize(1.5);
         hCOP->SetTitle(tit);
         hCOP->Draw("e1");
         hCOPbg->Draw("hist same");
         hCOPsig->Draw("hist same");
         hCOPtot->Draw("hist same");
         if (kUseCop)
         {
            lCoplo->Draw("c same");
            lCopup->Draw("c same");
         }

         c2->cd(3);
         hIM->GetXaxis()->SetRangeUser(fFitIMLo,fFitIMUp);
         hIM->GetYaxis()->SetRangeUser(rIM[0],rIM[1]);
         hIM->GetXaxis()->SetTitleSize(1.5);
         hIM->SetTitle(tit);
         hIM->Draw("e1");
         hIMbg->Draw("hist same");
         hIMsig->Draw("hist same");
         hIMtot->Draw("hist same");
         if (kUseIM)
         {
            lIMlo->Draw("c same");
            lIMup->Draw("c same");
         }

         c2->cd(4);
         legA->Draw();

         c2->cd(5);
         if (kUseMM) hMMsig->SetTitle(Form("int=%4.3e",intMM));
         else hMMsig->SetTitle("");
         hMMsig->GetXaxis()->SetTitleSize(1.5);
         hMMsig->GetXaxis()->SetRangeUser(fFitMMLo,fFitMMUp);
         hMMsig->GetYaxis()->SetRangeUser(r2MM[0],r2MM[1]);
         hMMsig->Draw("e1");
         hMMDiff->Draw("hist same");
         if (kUseMM)
         {
            lMM2lo->Draw("c same");
            lMM2up->Draw("c same");
         }

         c2->cd(6);
         if (kUseCop) hCOPsig->SetTitle(Form("int=%4.3e",intCop));
         else hCOPsig->SetTitle("");
         hCOPsig->GetXaxis()->SetTitleSize(1.5);
         hCOPsig->GetXaxis()->SetRangeUser(fFitCopLo,fFitCopUp);
         hCOPsig->GetYaxis()->SetRangeUser(r2Cop[0],r2Cop[1]);
         hCOPsig->Draw("e1");
         hCOPDiff->Draw("hist same");
         if (kUseCop)
         {
            lCop2lo->Draw("c same");
            lCop2up->Draw("c same");
         }

         c2->cd(7);
         if (kUseIM) hIMsig->SetTitle(Form("int=%4.3e",intIM));
         else hIMsig->SetTitle("");
         hIMsig->GetXaxis()->SetTitleSize(1.5);
         hIMsig->GetXaxis()->SetRangeUser(fFitIMLo,fFitIMUp);
         hIMsig->GetYaxis()->SetRangeUser(r2IM[0],r2IM[1]);
         hIMsig->Draw("e1");
         hIMDiff->Draw("hist same");
         if (kUseIM)
         {
            lIM2lo->Draw("c same");
            lIM2up->Draw("c same");
         }

         c2->cd(8);
         legB->Draw();

         fOut->cd();
         c2->Write();
      }
   }

   dfMM->GetData()->Scale(nCTBins/4./TMath::Pi());
   dfCop->GetData()->Scale(nCTBins/4./TMath::Pi());
   dfIM->GetData()->Scale(nCTBins/4./TMath::Pi());

   dfMM->GetData()->Scale(1.e6);
   dfCop->GetData()->Scale(1.e6);
   dfIM->GetData()->Scale(1.e6);

   TCanvas* c4 = new TCanvas("c4","c4",1300,1000);
   TGraphErrors* gfMM  = dfMM->CreateTotalIntegralSum();
   TGraphErrors* gfCop = dfCop->CreateTotalIntegralSum();
   TGraphErrors* gfIM  = dfIM->CreateTotalIntegralSum();
   gfMM->SetMarkerStyle(20);
   gfMM->SetMarkerColor(kBlue);
   gfMM->SetMarkerSize(1.5);
   gfCop->SetMarkerStyle(20);
   gfCop->SetMarkerColor(kGreen);
   gfMM->SetMarkerSize(1.5);
   gfIM->SetMarkerStyle(20);
   gfIM->SetMarkerColor(kRed);
   gfIM->SetMarkerSize(1.5);

   csl = new TLegend(0.7,0.7,0.85,0.85);
   csl->SetBorderSize(0);
   csl->SetFillColor(0);
   csl->SetTextSize(0.045);
   csl->SetTextFont(42);
   csl->AddEntry(gRef,"ref CS","p");
   csl->AddEntry(gfMM,"MM","p");
   csl->AddEntry(gfCop,"Cop","p");
   csl->AddEntry(gfIM,"IM","p");

   c4->cd();
   gRef->Draw("ap");
   gRef->GetYaxis()->SetRangeUser(-0.5,1.5*TMath::MaxElement(gRef->GetN(),gRef->GetY()));//1.5*gRef->GetMaximum());
   gfMM->Draw("psame");
   gfCop->Draw("psame");
   gfIM->Draw("psame");
   csl->Draw("same");

   c4->Write();

   fOut->cd();
   dBG->Write();
   dSig->Write();
   dfBG->Write();
   dfSig->Write();
   dMM->Write();
   dCop->Write();
   dIM->Write();
   dfMM->Write();
   dfCop->Write();
   dfIM->Write();
   fOut->Close();

   delete fOut;

   for (Int_t i = 0; i < nCTBins; i++)
   {
      delete hMM2d[i];
      delete hCop2d[i];
      delete hIM2d[i];
      delete hMMbg2d[i];
      delete hCopbg2d[i];
      delete hIMbg2d[i];
      delete hMMsig2d[i];
      delete hCopsig2d[i];
      delete hIMsig2d[i];
   }

   delete [] hMM2d;
   delete [] hCop2d;
   delete [] hIM2d;
   delete [] hMMbg2d;
   delete [] hCopbg2d;
   delete [] hIMbg2d;
   delete [] hMMsig2d;
   delete [] hCopsig2d;
   delete [] hIMsig2d;

   return;
}

//_______________________________________________________________________________
Double_t* TMCombinedFit::GetMinMaxRange(TH1* hA, TH1* hB)
{
   Double_t* m = new Double_t[2];

   Double_t minA = hA->GetBinContent(hA->GetMinimumBin());
   Double_t minB = hB->GetBinContent(hB->GetMinimumBin());
   Double_t maxA = hA->GetBinContent(hA->GetMaximumBin());
   Double_t maxB = hB->GetBinContent(hB->GetMaximumBin());

   m[0] = minA;
   m[1] = maxA;

   if (minA > minB)
      m[0] = minB;

   if (maxA < maxB)
      m[1] = maxB;

   m[0] = m[0] - 0.1*TMath::Abs(m[0]);
   m[1] = m[1] + 0.1*TMath::Abs(m[1]);

   return m;
}

//_______________________________________________________________________________
void TMCombinedFit::DivideFlux(TH2* h, TH1* Flux)
{
    for(Int_t i = 1; i <= h->GetNbinsX(); i++)
    {
        for(Int_t j = 1; j <= h->GetNbinsY(); j++)
        {
           Double_t Content = h->GetBinContent(i,j);
           Double_t FluxContent = Flux->GetBinContent(i);
           if(FluxContent>0.)
           {
               h->SetBinError(i, j, h->GetBinError(i,j)/FluxContent);
               h->SetBinContent(i, j, Content / FluxContent);
           }
           else
           {
               h->SetBinError(i, j, 0.);
               h->SetBinContent(i, j, 0.);
           }
        }
    }

    Printf("Divided %s with flux %s", h->GetName(), Flux->GetName());

    return;
}
//_______________________________________________________________________________
void TMCombinedFit::DivideEfficiency(TH2** h, TOEnergyThetaData* eff)
{
    // Divide the 2-dim signal histogram 'sig' by the efficiencies stored in 'eff'.
    // If 'nCTBins' is zero use class member 'fNCTBins' as number of cos(theta) bins,
    // otherwise use 'nCTBins'.

    // check number of energy bins
    if (h[0]->GetNbinsX() != eff->GetNEnergyBin())
    {
        Error("DivideEfficiency", "Numbers of energy bins in signal histogram and efficiency differ!");
        return;
    }

    // check number of cos(theta) bins
    if (nCTBins != eff->GetNCosThetaBin())
    {
        Error("DivideEfficiency", "Numbers of cos(theta) bins in signal histogram and efficiency differ!");
        return;
    }

    // loop over energy bins
    for (Int_t i = 0; i < h[0]->GetNbinsX(); i++)
    {
        // loop over cos(theta) bins
        for (Int_t j = 0; j < nCTBins; j++)
        {
            Double_t e = eff->GetValue(i, j);
            if (e < 1e-6) e = 1;

            // loop over y-axis
            for (Int_t k = 0; k < h[0]->GetNbinsY(); k++)
            {
                h[j]->SetBinContent(i+1, k+1, h[j]->GetBinContent(i+1, k+1) / e);
                h[j]->SetBinError(i+1, k+1, h[j]->GetBinError(i+1, k+1) / e);
            }
        }
    }

    Printf("Divided %s with efficiency %s",h[0]->GetName(), eff->GetName());

   return;
}
