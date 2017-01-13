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


#include "TMCuts.h"

//______________________________________________________________________________
void TMCuts::FitMMnew(TMObjectCollection* m, const Char_t* s, Int_t iCT, Int_t iRebinMM, TOEnergyThetaCut* e)
{
    Char_t szName[256];

    TFile* fData = 0;
    TFile* fMC = 0;
    TFile** fBG = 0;

    // check if mm fit func file for fit range already exists
//    if (!gSystem->AccessPathName("cuts/mm_fit.root"))
//    {
//       fcut = new TFile("cuts/mm_fit.root");
//       ecut = (TOEnergyThetaCut*)fcut->Get("mm_cut_func")->Clone("cut");
//       fcut->Close();
//    }

    // get number of BG files and init file array
    Int_t n = m->GetArraySize();
    fBG = new TFile*[n];

    // check if collection contains files
    TIter next(m->GetCollection());
    TObject* o;
    while ((o = (TObject*)next()))
    {
       if (TMTools::IsObjectType(o,"TFile*"))
       {
          Error("TMCuts::FitMM","Collection Objects must be of type TFile");
          gSystem->Exit(0);
       }
    }

    // get data object
    if (!(fData=(TFile*)m->GetDataObject()))
       Error("TMCuts::FitMM","No Data Object found");
    else
       Info("TMCuts::FitMM","Found Data Object %s",fData->GetName());

    // get mc signal object
    if (!(fMC=(TFile*)m->GetSignalObject()))
    {
       Error("TMCuts::FitMM","No MC Signal Object found");
    }
    else
       Info("TMCuts::FitMM","Found MC Signal Object %s",fMC->GetName());

    // get mc BG object
    if (!(fBG[0]=(TFile*)m->GetBackgroundObject()) && !(fBG=(TFile**)m->GetArrayObject()))
       Error("TMCuts::FitMM","No MC Background Objects found");
    else
       for (Int_t i = 0; i < n; i++)
          Info("TMCuts::FitMM","Found MC Background Object %s",fBG[i]->GetName());

    // prompt CT bin number
    Info("TMCuts::FitMM","Using %i CT Bins",iCT);

    // init cut object
    TOEnergyThetaCut* cut = new TOEnergyThetaCut(iCT, "mm_cut");

    // loop over CT bins
    for (Int_t i = 0; i < iCT; i++)
    {
       // get lower, upper cut funcs if available
       TF1* fLower = 0;
       TF1* fUpper = 0;
       if (e)
       {
          fLower = e->GetLowerFunction(i);
          fUpper = e->GetUpperFunction(i);
       }

       // check if histograms exist in files
       sprintf(szName,"%s_%i",s,i);
       if (!TMTools::CheckObject(szName,fData))
          Error("TMCuts::FitMM","Histogram %s in file %s not found", szName,fData->GetName());
       else
          Info("TMCuts::FitMM","Found Histogram %s in file %s", szName,fData->GetName());

       if (fMC)
       {
          if (!TMTools::CheckObject(szName,fMC))
             Error("TMCuts::FitMM","Histogram %s in file %s not found", szName,fMC->GetName());
          else
             Info("TMCuts::FitMM","Found Histogram %s in file %s", szName,fMC->GetName());
       }
       else
          Error("TMCuts::FitMM","No Signal Histograms found");
       
       for (Int_t j = 0; j < n; j++)
       {
          if (!TMTools::CheckObject(szName,fBG[j]))
             Error("TMCuts::FitMM","Histogram %s in file %s not found", szName,fBG[j]->GetName());
          else
             Info("TMCuts::FitMM","Found Histogram %s in file %s", szName,fBG[j]->GetName());
       }

       // get histograms
       TH2F* hData = (TH2F*)fData->Get(szName)->Clone("hData"); 
       TH2F* hMC = (TH2F*)fMC->Get(szName)->Clone("hMC"); 
       TH2F** hBG;
       hBG = new TH2F*[n];
       for (Int_t j = 0; j < n; j++)
          hBG[j] = (TH2F*)fBG[j]->Get(szName)->Clone(Form("hBG_%i",j));

       // init mean and sigma graphs
       TGraphErrors* gMean  = new TGraphErrors(hData->GetNbinsX());
       TGraphErrors* gSigma = new TGraphErrors(hData->GetNbinsX());

       Double_t pTmp[n+1];
       for (Int_t j = 0; j < n+1; j++)
          pTmp[j] = 1.;

       // loop over energy bins
       for (Int_t j = 0; j < hData->GetNbinsX(); j++)
       {
          // get MM projection histo
          TH1D* hDataP = (TH1D*)hData->ProjectionY(Form("hDataP_%i_%i",i+1,j+1),j+1,j+1,"e");
          hDataP->Rebin(iRebinMM);

          // init fitting histo array
          TH1** hFit;
          hFit = new TH1*[n+1];

          // pass mc signal histo
          hFit[0] = (TH1D*)hMC->ProjectionY(Form("hMCP_%i_%i",i+1,j+1),j+1,j+1,"e");
          hFit[0]->Rebin(iRebinMM);
          // pass mc bg histos
          for (Int_t k = 0; k < n; k++)
          {
             hFit[k+1] = (TH1D*)hBG[k]->ProjectionY(Form("hBGP_%i_%i_%i",k,i+1,j+1),j+1,j+1,"e");
             hFit[k+1]->Rebin(iRebinMM);
          }

          // fit with mc signal only first   
          TMFunctor::InitFunctor(1,hFit);
          TF1* fFitTmp = new TF1("fFitTmp",TMFunctor::MyFitFuncMC,hDataP->GetXaxis()->GetBinLowEdge(1),hDataP->GetXaxis()->GetBinUpEdge(hDataP->GetNbinsX()),1);
          fFitTmp->SetParameter(0,hDataP->GetMaximum()/hFit[0]->GetMaximum());
//          fFitTmp->SetParLimits(0,0.9*hDataP->GetMaximum()/hFit[0]->GetMaximum(),1.1*hDataP->GetMaximum()/hFit[0]->GetMaximum());

          // fit with mm cut positions as fit range if available
          if (e)
          {
             hDataP->Fit(fFitTmp,"+RW0MQ","",fLower->Eval(hDataP->GetXaxis()->GetBinCenter(j+1)),
                                             fUpper->Eval(hDataP->GetXaxis()->GetBinCenter(j+1)));
             Info("FitMM","Fitting Signal from %5.3f to %5.3f",
                                     fLower->Eval(hDataP->GetXaxis()->GetBinCenter(j+1)),
                                     fUpper->Eval(hDataP->GetXaxis()->GetBinCenter(j+1)));
          }
          else
          {
             hDataP->Fit(fFitTmp,"+RW0MQ","",-1000,1000);
             Info("FitMM","Fitting Signal in standard range");
          }

          // scale signal object
          Double_t p0 = fFitTmp->GetParameter(0);
          hFit[0]->Scale(p0);

          // fit also with mc BG now
          TF1* fFit = 0;
          TMFunctor::InitFunctor(n+1,hFit);
          fFit = new TF1("fFit",TMFunctor::MyFitFuncMC,hDataP->GetXaxis()->GetBinLowEdge(1),hDataP->GetXaxis()->GetBinUpEdge(hDataP->GetNbinsX()),n+1);
          Double_t pFit[n+1];
          pFit[0] = p0;
          for (Int_t k = 1; k < n+1; k++)
             pFit[k] = pTmp[k];
          fFit->SetParameters(pFit);
          fFit->SetParLimits(0,0.95*p0,1.05*p0);
          // accept only signal variations between 80-120%
//          fFit->SetParLimits(0,0.000001,100000.);
//          for (Int_t k = 0; k < n+1; k++)
//          {
//             if (j > 0)
//                fFit->SetParLimits(k,0.5*pTmp[k],1.5*pTmp[k]);
//             if (k==0)
//                fFit->SetParLimits(k,0.8*pFit[k],1.2*pFit[k]);
//                fFit->SetParLimits(k,1.e-3,1.e5);
//             else
//                fFit->SetParLimits(k,1.e-3,1.e5);
//          }

          hDataP->Fit(fFit,"+RW0MQ");

          for (Int_t k = 0; k < n+1; k++)
             pTmp[k] = fFit->GetParameter(k);

          // scale mc histos with fit parameters and add up all histos to total histo
          hFit[0]->Scale(fFit->GetParameter(0));
          TH1D* hTot = (TH1D*)hFit[0]->Clone(Form("hTot_%i_%i",i+1,j+1));
          for (Int_t k = 1; k < n+1; k++)
          {
             hFit[k]->Scale(fFit->GetParameter(k));
             hTot->Add(hFit[k]);
          }
          // add up BG histos to total BG histo
          TH1D* hTotBG = (TH1D*)hFit[1]->Clone(Form("hTotBG_%i_%i",i+1,j+1));
          for (Int_t k = 2; k < n+1; k++)
             hTotBG->Add(hFit[k]);

          // fit mc signal with gaussian to get lower, upper cut positions
          TF1* fGauss = new TF1(Form("fGauss_%i_%i",i+1,j+1), "gaus", -100, 100);
          fGauss->SetTitle(Form("%5.3f",hDataP->GetXaxis()->GetBinCenter(j+1)));
          fGauss->SetRange(-500.,500.);
//          fGauss->SetRange(TMTools::GetFirstNotEmptyBin(hFit[0]),TMTools::GetLastNotEmptyBin(hFit[0]));
          fGauss->SetParameters(1000.,0.,50.);
          fGauss->SetParLimits(0,0.,100000.);
	
          fGauss->SetParLimits(1,0.98*hFit[0]->GetBinCenter(hFit[0]->GetMaximumBin()),1.02*hFit[0]->GetBinCenter(hFit[0]->GetMaximumBin()));
          fGauss->SetParLimits(2,10.,100000.);
          hFit[0]->Fit(fGauss,"+RB0MQ");

          // add objects to collection
          TMObjectCollection* mm = new TMObjectCollection(Form("mmObj_%i_%i",i+1,j+1),Form("E_{#gamma}=%i MeV",(Int_t)hData->GetXaxis()->GetBinCenter(j+1)));
          mm->AddDataObject(hDataP);
//          TH1* hSig = (TH1*)hFit[0]->Clone(Form("hSig_%i_%i",i+1,j+1)); 
          mm->AddSignalObject(hFit[0]);
          mm->AddArrayObject(n,(TObject**)&hFit[1]);
          mm->AddBackgroundObject(hTotBG);
          mm->AddTotalObject(hTot);
          mm->AddFitObject(fGauss);
          mm->Write();

          // add mean and sigma cut positions
          gMean->SetPoint(j,hData->GetXaxis()->GetBinCenter(j+1),fGauss->GetParameter(1));
          gSigma->SetPoint(j,hData->GetXaxis()->GetBinCenter(j+1),fGauss->GetParameter(2));
          gMean->SetPointError(j,0.,0.);
          gSigma->SetPointError(j,0.,0.);

          delete [] hFit;
          delete fFit;
          delete fGauss;
          delete hDataP;
          delete hTot;
       }

       delete [] hBG;
       delete hData;
       delete hMC;

       // write mean and sigma cut positions
       gMean->Write(Form("Mean_%d", i));
       cut->SetMeanGraph(i, gMean);
       gSigma->Write(Form("Sigma_%d", i));
       cut->SetSigmaGraph(i, gSigma);

       delete gMean;
       delete gSigma;

    }
    delete [] fBG;

    // write cut positions
    cut->Write();
    delete cut;
}
//______________________________________________________________________________
void TMCuts::FitMM(TMObjectCollection* m, const Char_t* s, Int_t iCT, Int_t iRebinMM, Bool_t kFit, Bool_t kEmpty, Bool_t k2Empty)
{
    Char_t szName[256];
    Char_t szzName[256];

    TFile* fData = 0;
    TFile* fEmpty = 0;
    TFile* fMC = 0;
    TFile* fMC1 = 0;
    TFile* fMC2 = 0;
    TFile** fBG = 0;
    Int_t n = m->GetArraySize();

    fBG = new TFile*[n];

    TIter next(m->GetCollection());
    TObject* o;
    while ((o = (TObject*)next()))
    {
       if (TMTools::IsObjectType(o,"TFile*"))
       {
          Error("TMCuts::FitMM","Collection Objects must be of type TFile");
          gSystem->Exit(0);
       }
    }

    if (!(fData=(TFile*)m->GetDataObject()))
       Error("TMCuts::FitMM","No Data Object found");
    else
       Info("TMCuts::FitMM","Found Data Object %s",fData->GetName());

    if (!(fEmpty=(TFile*)m->GetEmptyTargetObject()))
       Error("TMCuts::FitMM","No Empty Target Object found");
    else
       Info("TMCuts::FitMM","Found Empty Target Object %s",fEmpty->GetName());

    if (!(fMC=(TFile*)m->GetSignalObject()))
    {
       if (!(fMC1=(TFile*)m->GetSignal1Object()))
          Error("TMCuts::FitMM","No MC Signal1 Object found");
       else
          Info("TMCuts::FitMM","Found MC Signal1 Object %s",fMC1->GetName());
       if (!(fMC2=(TFile*)m->GetSignal2Object()))
          Error("TMCuts::FitMM","No MC Signal2 Object found");
       else
          Info("TMCuts::FitMM","Found MC Signal2 Object %s",fMC2->GetName());
    }
    else
       Info("TMCuts::FitMM","Found MC Signal Object %s",fMC->GetName());
    if (!(fBG[0]=(TFile*)m->GetBackgroundObject()) && !(fBG=(TFile**)m->GetArrayObject()))
       Error("TMCuts::FitMM","No MC Background Objects found");
    else
       for (Int_t i = 0; i < n; i++)
          Info("TMCuts::FitMM","Found MC Background Object %s",fBG[i]->GetName());

    Info("TMCuts::FitMM","Using %i CT Bins",iCT);

    TOEnergyThetaCut* cut = 0;

    if (kFit) cut = new TOEnergyThetaCut(iCT, "mm_cut");

    for (Int_t i = 0; i < iCT; i++)
    {
       sprintf(szzName,"%s_%i",s,i);
       if (!TMTools::CheckObject(szzName,fData))
          Error("TMCuts::FitMM","Histogram %s in file %s not found", szzName,fData->GetName());
       else
          Info("TMCuts::FitMM","Found Histogram %s in file %s", szzName,fData->GetName());

       if (!TMTools::CheckObject(szzName,fEmpty))
          Error("TMCuts::FitMM","Histogram %s in file %s not found", szzName,fEmpty->GetName());
       else
          Info("TMCuts::FitMM","Found Histogram %s in file %s", szzName,fEmpty->GetName());

       if (fMC)
       {
          if (!TMTools::CheckObject(szzName,fMC))
             Error("TMCuts::FitMM","Histogram %s in file %s not found", szzName,fMC->GetName());
          else
             Info("TMCuts::FitMM","Found Histogram %s in file %s", szzName,fMC->GetName());
       }
       else if (fMC1 && fMC2)
       {
          if (!TMTools::CheckObject(szzName,fMC1))
             Error("TMCuts::FitMM","Histogram %s in file %s not found", szzName,fMC1->GetName());
          else
             Info("TMCuts::FitMM","Found Histogram %s in file %s", szzName,fMC1->GetName());
          if (!TMTools::CheckObject(szzName,fMC2))
             Error("TMCuts::FitMM","Histogram %s in file %s not found", szzName,fMC2->GetName());
          else
             Info("TMCuts::FitMM","Found Histogram %s in file %s", szzName,fMC2->GetName());
       }
       else
          Error("TMCuts::FitMM","No Signal Histograms found");
       for (Int_t j = 0; j < n; j++)
       {
          if (!TMTools::CheckObject(szzName,fBG[j]))
             Error("TMCuts::FitMM","Histogram %s in file %s not found", szzName,fBG[j]->GetName());
          else
             Info("TMCuts::FitMM","Found Histogram %s in file %s", szzName,fBG[j]->GetName());
       }

       TH2F* hData = 0;
       TH2F* hEmp = 0;

       if (kEmpty)
       {
          hData = SubtractEmptyTarget(fData, fEmpty, szzName);
       }
       else if (k2Empty)
       {
          hData = (TH2F*)fData->Get(szzName)->Clone("hData"); 
          hEmp  = (TH2F*)fEmpty->Get(szzName)->Clone("hEmp");
       }
       else
          hData = (TH2F*)fData->Get(szzName)->Clone("hData"); 
       
       TH2F* hMC = 0;
       TH2F* hMC1 = 0;
       TH2F* hMC2 = 0;

       if (fMC)
          hMC = (TH2F*)fMC->Get(szzName)->Clone("hMC"); 
       else 
       {
          hMC1 = (TH2F*)fMC1->Get(szzName)->Clone("hMC1"); 
          hMC2 = (TH2F*)fMC2->Get(szzName)->Clone("hMC2"); 
       }

       TH2F** hBG;
       hBG = new TH2F*[n];
    
       for (Int_t j = 0; j < n; j++)
       {
          sprintf(szName,"hBG_%i",j);
          hBG[j] = (TH2F*)fBG[j]->Get(szzName)->Clone(szName);
       }

       TGraphErrors* gMean = 0;
       TGraphErrors* gSigma = 0;

       TGraphErrors** gPar = 0;

       if (k2Empty)
       {
          gPar = new TGraphErrors*[3];

          for (Int_t j = 0; j < 3; j++)
             gPar[j] = new TGraphErrors(hData->GetNbinsX());
       }

       if (kFit)
       {
          gMean  = new TGraphErrors(hData->GetNbinsX());
          gSigma = new TGraphErrors(hData->GetNbinsX());
       }

       for (Int_t j = 0; j < hData->GetNbinsX(); j++)
       {
          sprintf(szName,"hDataP_%i_%i",i+1,j+1);
          TH1D* hDataP = (TH1D*)hData->ProjectionY(szName,j+1,j+1,"e");
          hDataP->Rebin(iRebinMM);

          TH1** hFit;
          if (hMC)
             hFit = new TH1*[n+1];
          else
             hFit = new TH1*[n+2];

          TF1* fFitTmp = 0;

          if (hMC)
          {
             sprintf(szName,"hMCP_%i_%i",i+1,j+1);
             hFit[0] = (TH1D*)hMC->ProjectionY(szName,j+1,j+1,"e");
             hFit[0]->Rebin(iRebinMM);
   
             TMFunctor::InitFunctor(1,hFit);
             fFitTmp = new TF1("fFitTmp",TMFunctor::MyFitFuncMC,hDataP->GetXaxis()->GetBinLowEdge(1),hDataP->GetXaxis()->GetBinUpEdge(hDataP->GetNbinsX()),1);
             fFitTmp->SetParameter(0,1.);
             fFitTmp->SetParLimits(0,0.001,100000.);
          }
          else
          {
             sprintf(szName,"hMC1P_%i_%i",i+1,j+1);
             hFit[0] = (TH1D*)hMC1->ProjectionY(szName,j+1,j+1,"e");
             hFit[0]->Rebin(iRebinMM);
             sprintf(szName,"hMC2P_%i_%i",i+1,j+1);
             hFit[1] = (TH1D*)hMC2->ProjectionY(szName,j+1,j+1,"e");
             hFit[1]->Rebin(iRebinMM);

             TMFunctor::InitFunctor(2,hFit);
             fFitTmp = new TF1("fFitTmp",TMFunctor::MyFitFuncMC,hDataP->GetXaxis()->GetBinLowEdge(1),hDataP->GetXaxis()->GetBinUpEdge(hDataP->GetNbinsX()),2);

             Double_t pFitTmp[2];
             for (Int_t k = 0; k < 2; k++)
                pFitTmp[k]=1.;
             fFitTmp->SetParameters(pFitTmp);

             for (Int_t k = 0; k < 2; k++)
             {
                fFitTmp->SetParLimits(k,0.001,100000.);
             }
          }

          hDataP->Fit(fFitTmp,"+RW0MQ");

          TF1* fFit = 0;
          Int_t nMax = 0;

          if (hMC)
          {
             sprintf(szName,"hMCP_%i_%i",i+1,j+1);
             hFit[0] = (TH1D*)hMC->ProjectionY(szName,j+1,j+1,"e");
             hFit[0]->Rebin(iRebinMM);
             for (Int_t k = 0; k < n; k++)
             {
                sprintf(szName,"hBGP_%i_%i_%i",k,i+1,j+1);
                hFit[k+1] = (TH1D*)hBG[k]->ProjectionY(szName,j+1,j+1,"e");
                hFit[k+1]->Rebin(iRebinMM);
             }
   
             TMFunctor::InitFunctor(n+1,hFit);
             fFit = new TF1("fFit",TMFunctor::MyFitFuncMC,hDataP->GetXaxis()->GetBinLowEdge(1),hDataP->GetXaxis()->GetBinUpEdge(hDataP->GetNbinsX()),n+1);
             Double_t pFit[n+1];
             pFit[0] = fFitTmp->GetParameter(0);
             for (Int_t k = 1; k < n+1; k++)
                pFit[k] = 1.;
             fFit->SetParameters(pFit);
             fFit->SetParLimits(0,0.1*pFit[0],1.5*pFit[0]);
             for (Int_t k = 1; k < n+1; k++)
             {
                fFit->SetParLimits(k,0.,100000.);
             }
             nMax = n+1;
          }
          else
          {
             sprintf(szName,"hMC1P_%i_%i",i+1,j+1);
             hFit[0] = (TH1D*)hMC1->ProjectionY(szName,j+1,j+1,"e");
             hFit[0]->Rebin(iRebinMM);
             sprintf(szName,"hMC2P_%i_%i",i+1,j+1);
             hFit[1] = (TH1D*)hMC2->ProjectionY(szName,j+1,j+1,"e");
             hFit[1]->Rebin(iRebinMM);

             for (Int_t k = 0; k < n; k++)
             {
                sprintf(szName,"hBGP_%i_%i_%i",k,i+1,j+1);
                hFit[k+2] = (TH1D*)hBG[k]->ProjectionY(szName,j+1,j+1,"e");
                hFit[k+2]->Rebin(iRebinMM);
             }

             TMFunctor::InitFunctor(n+2,hFit);
             fFit = new TF1("fFit",TMFunctor::MyFitFuncMC,hDataP->GetXaxis()->GetBinLowEdge(1),hDataP->GetXaxis()->GetBinUpEdge(hDataP->GetNbinsX()),n+2);
             Double_t pFit[n+2];
             pFit[0] = fFitTmp->GetParameter(0);
             pFit[1] = fFitTmp->GetParameter(1);
             for (Int_t k = 2; k < n+2; k++)
                pFit[k]=10000.;
             fFit->SetParameters(pFit);

             fFit->SetParLimits(0,0.8*pFit[0],1.2*pFit[0]);
             fFit->SetParLimits(1,0.8*pFit[1],1.2*pFit[1]);
             for (Int_t k = 2; k < n+2; k++)
             {
                fFit->SetParLimits(k,0.,100000.);
             }
             nMax = n+2;
          }

          hDataP->Fit(fFit,"+RW0MQ");

          hFit[0]->Scale(fFit->GetParameter(0));
          sprintf(szName,"hTot_%i_%i",i+1,j+1);
          TH1D* hTot = (TH1D*)hFit[0]->Clone(szName);

//printf("Height Bin 80: %s: %5.3f  factor: %5.3f\n", hFit[0]->GetName(),hFit[0]->GetBinContent(72),fFit->GetParameter(0));
//printf("Height Bin 80: %s: %5.3f\n", hTot->GetName(),hTot->GetBinContent(72));
          for (Int_t k = 1; k < nMax; k++)
          {
             hFit[k]->Scale(fFit->GetParameter(k));
//printf("Height Bin 80: %s: %5.3f  factor: %5.3f\n", hFit[k]->GetName(),hFit[k]->GetBinContent(72),fFit->GetParameter(k));
             hTot->Add(hFit[k]);
//printf("Height Bin 80: %s: %5.3f\n", hTot->GetName(),hTot->GetBinContent(72));

          }

          sprintf(szName,"mmObj_%i_%i",i+1,j+1);
          sprintf(szzName,"E_{#gamma}=%i MeV",(Int_t)hData->GetXaxis()->GetBinCenter(j+1));
          TMObjectCollection* mm = new TMObjectCollection(szName,szzName);
          mm->AddDataObject(hDataP);

          TH1* hSig = 0;
          if (hMC)
          {
             sprintf(szName,"hSig_%i_%i",i+1,j+1);
             hSig = (TH1*)hFit[0]->Clone(szName); 
             mm->AddSignalObject(hSig);
             mm->AddArrayObject(n,(TObject**)&hFit[1]);
          }
          else
          {
             sprintf(szName,"hSig_%i_%i",i+1,j+1);
             hSig = (TH1*)hFit[0]->Clone(szName);
             hSig->Add(hFit[1]);
             mm->AddSignalObject(hSig);
             mm->AddArrayObject(n,(TObject**)&hFit[2]);
          }
          mm->AddTotalObject(hTot);

          TH1D* hTotBG = 0;
          sprintf(szName,"hTotBG_%i_%i",i+1,j+1);

          if (hMC)
          {
             hTotBG = (TH1D*)hFit[1]->Clone(szName);

             for (Int_t k = 2; k < nMax; k++)
                hTotBG->Add(hFit[k]);
          }
          else
          {
             hTotBG = (TH1D*)hFit[2]->Clone(szName);

             for (Int_t k = 3; k < nMax; k++)
                hTotBG->Add(hFit[k]);
          }

          TF1* fGauss = 0;
          if (kFit)
          {
             sprintf(szName,"fGauss_%i_%i",i+1,j+1);
             fGauss = new TF1(szName,
                                   "gaus",
                                   -1000,
                                   1000);
             sprintf(szName,"%5.3f",hDataP->GetXaxis()->GetBinCenter(j+1));
             fGauss->SetTitle(szName);
             fGauss->SetRange(TMTools::GetFirstNotEmptyBin(hFit[0]),TMTools::GetLastNotEmptyBin(hFit[0]));
             fGauss->SetParameters(1000.,0.,50.);
             fGauss->SetParLimits(0,0.,100000.);
             fGauss->SetParLimits(2,10.,100000.);
             if (hMC)
                hFit[0]->Fit(fGauss,"+RB0MQ");
             else
                hSig->Fit(fGauss,"+RB0MQ");
             mm->AddFitObject(fGauss);
             gMean->SetPoint(j,hData->GetXaxis()->GetBinCenter(j+1),fGauss->GetParameter(1));
             gSigma->SetPoint(j,hData->GetXaxis()->GetBinCenter(j+1),fGauss->GetParameter(2));
//             gMean->SetPointError(j,hData->GetXaxis()->GetBinWidth(i+1)/2.,fGauss->GetParError(1));
//             gSigma->SetPointError(j,hData->GetXaxis()->GetBinWidth(i+1)/2.,fGauss->GetParError(2));
             gMean->SetPointError(j,0.,0.);
             gSigma->SetPointError(j,0.,0.);
          }

          mm->Write();

//-----------------------------------------------
// second fit
          if (k2Empty)
          {
             sprintf(szName,"hEmpP_0_%i_%i",i+1,j+1);
             TH1D* hEmptyP = (TH1D*)hEmp->ProjectionY(szName,j+1,j+1,"e");
             hEmptyP->Rebin(iRebinMM);

             TH1** h2Fit;
             h2Fit = new TH1*[3];
  
//             h2Fit[0] = hTot;
             h2Fit[0] = hSig;
             h2Fit[1] = hTotBG;
             h2Fit[2] = hEmptyP;
   
             TMFunctor::InitFunctor(3,h2Fit);
             TF1* f2Fit = new TF1("f2Fit",TMFunctor::MyFitFuncEmpty,hDataP->GetXaxis()->GetBinLowEdge(1),hDataP->GetXaxis()->GetBinUpEdge(hDataP->GetNbinsX()),3);
             Double_t p2Fit[3];
             for (Int_t k = 0; k < 3; k++)
                p2Fit[k]=1.;
             f2Fit->SetParameters(p2Fit);
             for (Int_t k = 0; k < 3; k++)
             {
                f2Fit->SetParLimits(k,0.,100000.);
             }
//             f2Fit->SetParLimits(0,0.9,1.1);
//             f2Fit->SetParLimits(1,0.9,1.1);
//             f2Fit->SetParLimits(2,0.,100000.);

             hDataP->Fit(f2Fit,"+RW0MQ");

             h2Fit[0]->Scale(f2Fit->GetParameter(0));
             h2Fit[1]->Scale(f2Fit->GetParameter(1));
             h2Fit[2]->Scale(f2Fit->GetParameter(2));

             gPar[0]->SetPoint(j,hData->GetXaxis()->GetBinCenter(j+1),f2Fit->GetParameter(0));
             gPar[0]->SetPointError(j,0.,f2Fit->GetParError(0));
             gPar[1]->SetPoint(j,hData->GetXaxis()->GetBinCenter(j+1),f2Fit->GetParameter(1));
             gPar[1]->SetPointError(j,0.,f2Fit->GetParError(1));
             gPar[2]->SetPoint(j,hData->GetXaxis()->GetBinCenter(j+1),f2Fit->GetParameter(2));
             gPar[2]->SetPointError(j,0.,f2Fit->GetParError(2));

             sprintf(szName,"h2Tot_%i_%i",i+1,j+1);
             TH1D* h2Tot = (TH1D*)h2Fit[0]->Clone(szName);

             h2Tot->Add(h2Fit[1]);
             h2Tot->Add(h2Fit[2]);
  
             sprintf(szName,"mm2Obj_%i_%i",i+1,j+1);
             sprintf(szzName,"E_{#gamma}=%i MeV",(Int_t)hData->GetXaxis()->GetBinCenter(j+1));
             TMObjectCollection* m2m = new TMObjectCollection(szName,szzName);
             m2m->AddDataObject(hDataP);

             m2m->AddSignalObject(h2Fit[0]);
             m2m->AddBackgroundObject(h2Fit[1]);
             m2m->AddArrayObject(1,(TObject**)&h2Fit[2]);
             m2m->AddTotalObject(h2Tot);

             m2m->Write();

             delete [] h2Fit;
             delete f2Fit;
             delete h2Tot;
          }
//-----------------------------------------------

//          if (hMC)
//          {
//             for (Int_t k = 0; k < n+1; k++)
//             {
//                printf("Here16_%i\n",k);
//                delete hFit[i];
//             }
//          }
//          else
//          {
//             for (Int_t k = 0; k < n+2; k++)
//             {
//                printf("Here16_%i\n",k);
//                delete hFit[i];
//             }
//          }
//          printf("Here16\n");

          delete [] hFit;
          delete fFit;
          if (kFit) delete fGauss;
          delete hDataP;
          delete hTot;
       }

       delete [] hBG;

       delete hData;

       if (hMC)
          delete hMC;
       else
       {
          delete hMC1;
          delete hMC2;
       }

       if (kFit)
       {  
//          TOHUtils::RebinGraph(gMean, 2);
//          TOHUtils::SmoothGraph(gMean);
          sprintf(szName, "Mean_%d", i);
          gMean->Write(szName);
          cut->SetMeanGraph(i, gMean);
//          TOHUtils::RebinGraph(gSigma, 2);
//          TOHUtils::SmoothGraph(gSigma);
          sprintf(szName, "Sigma_%d", i);
          gSigma->Write(szName);
          cut->SetSigmaGraph(i, gSigma);

          delete gMean;
          delete gSigma;
       }

       if (k2Empty)
       {
          sprintf(szName,"FitParameter_0_%i",i);
          gPar[0]->Write(szName);
          sprintf(szName,"FitParameter_1_%i",i);
          gPar[1]->Write(szName);
          sprintf(szName,"FitParameter_2_%i",i);
          gPar[2]->Write(szName);

//          delete gPar[0];
//          delete gPar[1];
//          delete gPar[2];
          delete hEmp;
       }
//       delete [] gPar;
    }
    delete [] fBG;

    if (kFit)
    {
       cut->Write();
       delete cut;
    }
}

//______________________________________________________________________________
void TMCuts::SubtractIMfromMM(Int_t iCT, TFile* fDataIn, TFile* fMCIn, const Char_t* szHin, const Char_t* szHout)
{
   Char_t szName[256];
   Char_t szzName[256];

//   sprintf(szName,"%s_subtracted",fDataIn->GetName());
//   TFile* f = new TFile(szName,"recreate");

//   f->cd();

   // loop over all CT histos
   for (Int_t i = 0; i < iCT; i++)
//   for (Int_t i = 4; i < 5; i++)
   {
      sprintf(szName,"%s_%i",szHin,i);
      sprintf(szzName,"hDataTmp_%i",i);
      TH3F* hDataTmp = (TH3F*)fDataIn->Get(szName)->Clone(szzName);
      sprintf(szzName,"hMCTmp_%i",i);
      TH3F* hMCTmp = (TH3F*)fMCIn->Get(szName)->Clone(szzName);

      sprintf(szName,"%s_%i_old",szHout,i);
      TH2D* hOld = (TH2D*)hDataTmp->Project3D("yx")->Clone(szName);
      hOld->Write();
      sprintf(szName,"%s_%i",szHout,i);
      TH2D* hNew = (TH2D*)hDataTmp->Project3D("yx")->Clone(szName);

      // loop over all TCbins (x-axis) of each histo
      for (Int_t j = 0; j < hDataTmp->GetNbinsX(); j++)
//      for (Int_t j = 25; j < 26; j++)
      {
         // loop over all MMbins (y-axis) of each histo
         for (Int_t k = 0; k < hDataTmp->GetNbinsY(); k++)
         {
            sprintf(szName,"hDataTmpP_%i_%i_%i",i+1,j+1,k+1);
            TH1D* hDataTmpP = (TH1D*)hDataTmp->ProjectionZ(szName,j+1,j+1,k+1,k+1);
            sprintf(szName,"hMCTmpP_%i_%i_%i",i+1,j+1,k+1);
            TH1D* hMCTmpP = (TH1D*)hMCTmp->ProjectionZ(szName,j+1,j+1,k+1,k+1);

            TH1** hFit;
            hFit = new TH1*[1];

            hFit[0] = hMCTmpP; 

            TMFunctor::InitFunctor(1,hFit,1,110.,160.);

            sprintf(szName,"fRej_%i_%i_%i",i+1,j+1,k+1);
            TF1* fRej = new TF1(szName,TMFunctor::Pol1Exc,0.,300.,2);
            fRej->SetParameters(2,-1);
            hDataTmpP->Fit(fRej,"0Q");

            TMFunctor::SetParameters(fRej->GetParameters());

            sprintf(szName,"fFitTmp_%i_%i_%i",i+1,j+1,k+1);
            TF1* fFit = new TF1(szName,TMFunctor::FitTH1Pol1,110.,160.,1);
            fFit->SetParameter(0,1.);
            hDataTmpP->Fit(fFit,"+RW0MQ");
            hMCTmpP->Scale(fFit->GetParameter(0));

            sprintf(szName,"fPol_%i_%i_%i",i+1,j+1,k+1);
            TF1* fPol = new TF1(szName,"pol1",0.,300.);
            fPol->SetParameters(fRej->GetParameters());
            fPol->SetLineColor(2);

            sprintf(szName,"hTot_%i_%i_%i",i+1,j+1,k+1);
            TH1D* hTot = (TH1D*)hMCTmpP->Clone(szName);
            hTot->Add(fPol);

            Double_t fCorrErrA = 0.;
            Double_t fCorrErrB = 0.;

            Double_t fCorrA = hMCTmpP->IntegralAndError(hMCTmpP->FindBin(110.),hMCTmpP->FindBin(160.),fCorrErrA);
            Double_t fCorrB = hDataTmpP->IntegralAndError(hDataTmpP->FindBin(110.),hMCTmpP->FindBin(160.),fCorrErrB);

            Double_t fCorr = 0.;
            Double_t fCorrErr = 0.;

            if (fCorrA > 0. && fCorrB > 0.)
            {
               fCorr    = fCorrA / fCorrB;
               fCorrErr = TMath::Sqrt(TMath::Power((fCorrErrA/fCorrB),2.)+TMath::Power(fCorrA*fCorrErrB/fCorrB/fCorrB,2.));
            }

            hNew->SetBinContent(j+1,k+1,hNew->GetBinContent(j+1,k+1)*fCorr);
            hNew->SetBinError(j+1,k+1,TMath::Sqrt(TMath::Power(fCorr*hNew->GetBinError(j+1,k+1),2.)+TMath::Power(hNew->GetBinContent(j+1,k+1)*fCorrErr,2.)));

            sprintf(szName,"imObj_%i_%i_%i",i+1,j+1,k+1);
            sprintf(szzName,"MM: %5.3f SB: %5.3f",hNew->GetYaxis()->GetBinCenter(k+1),fCorr);
            TMObjectCollection* oIM = new TMObjectCollection(szName,szzName);
            oIM->AddDataObject(hDataTmpP);
            oIM->AddSignalObject(hMCTmpP);
            oIM->AddBackgroundObject(fPol);
            oIM->AddTotalObject(hTot);
//
            oIM->Write();

            delete fFit;
            delete hFit[0];
            delete [] hFit;
//            delete oIM;
         }

         Info("TMCuts::SubtractIMfromMM","Finished CTbin %i and tCbin %i",i+1,j+1);
      }

      hNew->Write();
   }

//   delete f;

   return;
}

//______________________________________________________________________________
TH2F* TMCuts::SubtractEmptyTarget(TFile* f, TFile* fCorr, const Char_t* szName)
{
   Char_t szzName[256];

   TH2F* hNewData = (TH2F*)f->Get(szName)->Clone("hNewData");
   TH2F* hEmpty   = (TH2F*)fCorr->Get(szName)->Clone("hEmpty");
   TH1D* hFluxDataTmp  = (TH1D*)f->Get("FluxEg")->Clone("hFluxDataTmp");
   TH1D* hFluxEmptyTmp = (TH1D*)fCorr->Get("FluxEg")->Clone("hFluxEmptyTmp");

   TH2F* hFluxData = new TH2F("hFluxData","hFluxData",
                              hFluxDataTmp->GetNbinsX(),
                              hFluxDataTmp->GetXaxis()->GetBinLowEdge(1),
                              hFluxDataTmp->GetXaxis()->GetBinUpEdge(hFluxDataTmp->GetNbinsX()),
                              hNewData->GetNbinsY(),
                              hNewData->GetYaxis()->GetBinLowEdge(1),
                              hNewData->GetYaxis()->GetBinUpEdge(hNewData->GetNbinsY()));

   TH2F* hFluxEmpty = new TH2F("hFluxEmpty","hFluxEmpty",
                              hFluxEmptyTmp->GetNbinsX(),
                              hFluxEmptyTmp->GetXaxis()->GetBinLowEdge(1),
                              hFluxEmptyTmp->GetXaxis()->GetBinUpEdge(hFluxEmptyTmp->GetNbinsX()),
                              hEmpty->GetNbinsY(),
                              hEmpty->GetYaxis()->GetBinLowEdge(1),
                              hEmpty->GetYaxis()->GetBinUpEdge(hEmpty->GetNbinsY()));

   for (Int_t i = 0; i < hFluxDataTmp->GetNbinsX(); i++)
   {
      for (Int_t j = 0; j < hNewData->GetNbinsY(); j++)
      {
         hFluxData->SetBinContent(i+1,j+1,hFluxDataTmp->GetBinContent(i+1));
         hFluxEmpty->SetBinContent(i+1,j+1,hFluxEmptyTmp->GetBinContent(i+1));
      }
   }

   hNewData->Divide(hFluxData); 
   hEmpty->Divide(hFluxEmpty); 

   Int_t ii = 0;
   sscanf(szName,"%*[^_]_%*[^_]_%i",&ii);

   for (Int_t i = 0; i < hNewData->GetNbinsX(); i++)
   {
      TH1** hE;
      hE = new TH1*[2];

      sprintf(szzName,"empty_0_%i_%i",ii+1,i+1);
      hE[0] = hNewData->ProjectionY(szzName,i+1,i+1);
      sprintf(szzName,"empty_1_%i_%i",ii+1,i+1);
      hE[1] = hEmpty->ProjectionY(szzName,i+1,i+1);
      hE[1]->SetLineColor(2);

      sprintf(szzName,"emptyObj_%i_%i",ii+1,i+1);
      TMObjectCollection* oEmpty = new TMObjectCollection(szzName,szzName);
   
      oEmpty->AddArrayObject(2,(TObject**)&hE[0]);
      oEmpty->Write();
    
//      delete hE[0];
//      delete hE[1];
//      delete oEmpty;
   }

   hNewData->Add(hEmpty,-1);
   hNewData->Multiply(hFluxData);

   delete hFluxData;
   delete hFluxEmpty;

   return hNewData;
}

//______________________________________________________________________________
void TMCuts::FitCopNew(TMObjectCollection* m, const Char_t* s, Int_t iCT, Int_t iRebinMM, TOEnergyThetaCut* e)
{
    Char_t szName[256];

    TFile* fData = 0;
    TFile* fMC = 0;
    TFile** fBG = 0;

    // check if mm fit func file for fit range already exists
//    if (!gSystem->AccessPathName("cuts/mm_fit.root"))
//    {
//       fcut = new TFile("cuts/mm_fit.root");
//       ecut = (TOEnergyThetaCut*)fcut->Get("mm_cut_func")->Clone("cut");
//       fcut->Close();
//    }

    // get number of BG files and init file array
    Int_t n = m->GetArraySize();
    fBG = new TFile*[n];

    // check if collection contains files
    TIter next(m->GetCollection());
    TObject* o;
    while ((o = (TObject*)next()))
    {
       if (TMTools::IsObjectType(o,"TFile*"))
       {
          Error("TMCuts::FitCop","Collection Objects must be of type TFile");
          gSystem->Exit(0);
       }
    }

    // get data object
    if (!(fData=(TFile*)m->GetDataObject()))
       Error("TMCuts::FitCop","No Data Object found");
    else
       Info("TMCuts::FitCop","Found Data Object %s",fData->GetName());

    // get mc signal object
    if (!(fMC=(TFile*)m->GetSignalObject()))
    {
       Error("TMCuts::FitCop","No MC Signal Object found");
    }
    else
       Info("TMCuts::FitCop","Found MC Signal Object %s",fMC->GetName());

    // get mc BG object
    if (!(fBG[0]=(TFile*)m->GetBackgroundObject()) && !(fBG=(TFile**)m->GetArrayObject()))
       Error("TMCuts::FitCop","No MC Background Objects found");
    else
       for (Int_t i = 0; i < n; i++)
          Info("TMCuts::FitCop","Found MC Background Object %s",fBG[i]->GetName());

    // prompt CT bin number
    Info("TMCuts::FitCop","Using %i CT Bins",iCT);

    // init cut object
    TOEnergyThetaCut* cut = new TOEnergyThetaCut(iCT, "cop_cut");

    // loop over CT bins
    for (Int_t i = 0; i < iCT; i++)
    {
       // get lower, upper cut funcs if available
       TF1* fLower = 0;
       TF1* fUpper = 0;
       if (e)
       {
          fLower = e->GetLowerFunction(i);
          fUpper = e->GetUpperFunction(i);
       }

       // check if histograms exist in files
       sprintf(szName,"%s_%i",s,i);
       if (!TMTools::CheckObject(szName,fData))
          Error("TMCuts::FitCop","Histogram %s in file %s not found", szName,fData->GetName());
       else
          Info("TMCuts::FitCop","Found Histogram %s in file %s", szName,fData->GetName());

       if (fMC)
       {
          if (!TMTools::CheckObject(szName,fMC))
             Error("TMCuts::FitCop","Histogram %s in file %s not found", szName,fMC->GetName());
          else
             Info("TMCuts::FitCop","Found Histogram %s in file %s", szName,fMC->GetName());
       }
       else
          Error("TMCuts::FitCop","No Signal Histograms found");
       
       for (Int_t j = 0; j < n; j++)
       {
          if (!TMTools::CheckObject(szName,fBG[j]))
             Error("TMCuts::FitCop","Histogram %s in file %s not found", szName,fBG[j]->GetName());
          else
             Info("TMCuts::FitCop","Found Histogram %s in file %s", szName,fBG[j]->GetName());
       }

       // get histograms
       TH2F* hData = (TH2F*)fData->Get(szName)->Clone("hData"); 
       TH2F* hMC = (TH2F*)fMC->Get(szName)->Clone("hMC"); 
       TH2F** hBG;
       hBG = new TH2F*[n];
       for (Int_t j = 0; j < n; j++)
          hBG[j] = (TH2F*)fBG[j]->Get(szName)->Clone(Form("hBG_%i",j));

       // init mean and sigma graphs
       TGraphErrors* gMean  = new TGraphErrors(hData->GetNbinsX());
       TGraphErrors* gSigma = new TGraphErrors(hData->GetNbinsX());

       Double_t pTmp[n+1];
       for (Int_t j = 0; j < n+1; j++)
          pTmp[j] = 1.;

       // loop over energy bins
       for (Int_t j = 0; j < hData->GetNbinsX(); j++)
       {
          // get MM projection histo
          TH1D* hDataP = (TH1D*)hData->ProjectionY(Form("hDataP_%i_%i",i+1,j+1),j+1,j+1,"e");
          hDataP->Rebin(iRebinMM);

          // init fitting histo array
          TH1** hFit;
          hFit = new TH1*[n+1];

          // pass mc signal histo
          hFit[0] = (TH1D*)hMC->ProjectionY(Form("hMCP_%i_%i",i+1,j+1),j+1,j+1,"e");
          hFit[0]->Rebin(iRebinMM);
          // pass mc bg histos
          for (Int_t k = 0; k < n; k++)
          {
             hFit[k+1] = (TH1D*)hBG[k]->ProjectionY(Form("hBGP_%i_%i_%i",k,i+1,j+1),j+1,j+1,"e");
             hFit[k+1]->Rebin(iRebinMM);
          }

          // fit with mc signal only first   
          TMFunctor::InitFunctor(1,hFit);
          TF1* fFitTmp = new TF1("fFitTmp",TMFunctor::MyFitFuncMC,hDataP->GetXaxis()->GetBinLowEdge(1),hDataP->GetXaxis()->GetBinUpEdge(hDataP->GetNbinsX()),1);
          fFitTmp->SetParameter(0,hDataP->GetMaximum()/hFit[0]->GetMaximum());
//          fFitTmp->SetParLimits(0,0.9*hDataP->GetMaximum()/hFit[0]->GetMaximum(),1.1*hDataP->GetMaximum()/hFit[0]->GetMaximum());

          // fit with mm cut positions as fit range if available
          if (e)
          {
             hDataP->Fit(fFitTmp,"+RW0MQ","",fLower->Eval(hDataP->GetXaxis()->GetBinCenter(j+1)),
                                             fUpper->Eval(hDataP->GetXaxis()->GetBinCenter(j+1)));
             Info("FitCop","Fitting Signal from %5.3f to %5.3f",
                                     fLower->Eval(hDataP->GetXaxis()->GetBinCenter(j+1)),
                                     fUpper->Eval(hDataP->GetXaxis()->GetBinCenter(j+1)));
          }
          else
          {
             hDataP->Fit(fFitTmp,"+RW0MQ","",140,220);
             Info("FitCop","Fitting Signal in standard range");
          }

          // scale signal object
          Double_t p0 = fFitTmp->GetParameter(0);
          hFit[0]->Scale(p0);

          // fit also with mc BG now
          TF1* fFit = 0;
          TMFunctor::InitFunctor(n+1,hFit);
          fFit = new TF1("fFit",TMFunctor::MyFitFuncMC,hDataP->GetXaxis()->GetBinLowEdge(1),hDataP->GetXaxis()->GetBinUpEdge(hDataP->GetNbinsX()),n+1);
          Double_t pFit[n+1];
          pFit[0] = p0;
          for (Int_t k = 1; k < n+1; k++)
             pFit[k] = pTmp[k];
          fFit->SetParameters(pFit);
          fFit->SetParLimits(0,0.95*p0,1.05*p0);
          // accept only signal variations between 80-120%
//          fFit->SetParLimits(0,0.000001,100000.);
//          for (Int_t k = 0; k < n+1; k++)
//          {
//             if (j > 0)
//                fFit->SetParLimits(k,0.5*pTmp[k],1.5*pTmp[k]);
//             if (k==0)
                fFit->SetParLimits(0,0.95*pFit[0],1.15*pFit[0]);
//                fFit->SetParLimits(k,1.e-3,1.e5);
//             else
//                fFit->SetParLimits(k,1.e-3,1.e5);
//          }

          hDataP->Fit(fFit,"+RW0MQ");

          for (Int_t k = 0; k < n+1; k++)
             pTmp[k] = fFit->GetParameter(k);

          // scale mc histos with fit parameters and add up all histos to total histo
          hFit[0]->Scale(fFit->GetParameter(0));
          TH1D* hTot = (TH1D*)hFit[0]->Clone(Form("hTot_%i_%i",i+1,j+1));
          for (Int_t k = 1; k < n+1; k++)
          {
             hFit[k]->Scale(fFit->GetParameter(k));
             hTot->Add(hFit[k]);
          }
          // add up BG histos to total BG histo
          TH1D* hTotBG = (TH1D*)hFit[1]->Clone(Form("hTotBG_%i_%i",i+1,j+1));
          for (Int_t k = 2; k < n+1; k++)
             hTotBG->Add(hFit[k]);

          // fit mc signal with gaussian to get lower, upper cut positions
          TF1* fGauss = new TF1(Form("fGauss_%i_%i",i+1,j+1), "gaus", 0, 360);
          fGauss->SetTitle(Form("%5.3f",hDataP->GetXaxis()->GetBinCenter(j+1)));
	  fGauss->SetRange(100.,260.);
	  fGauss->SetParameters(hFit[0]->GetMaximum(),180.,10.);
	  fGauss->SetParLimits(0,0.8*hFit[0]->GetMaximum(),1.2*hFit[0]->GetMaximum());
	  fGauss->SetParLimits(1,170.,190.);
	  fGauss->SetParLimits(2,5.,50.);
//          fGauss->SetRange(TMTools::GetFirstNotEmptyBin(hFit[0]),TMTools::GetLastNotEmptyBin(hFit[0]));
//          fGauss->SetParameters(1000.,180.,50.);
//          fGauss->SetParLimits(0,0.,100000.);
//          fGauss->SetParLimits(2,10.,100000.);
          hFit[0]->Fit(fGauss,"+RB0MQ");

          // add objects to collection
          TMObjectCollection* mm = new TMObjectCollection(Form("copObj_%i_%i",i+1,j+1),Form("E_{#gamma}=%i MeV",(Int_t)hData->GetXaxis()->GetBinCenter(j+1)));
          mm->AddDataObject(hDataP);
//          TH1* hSig = (TH1*)hFit[0]->Clone(Form("hSig_%i_%i",i+1,j+1)); 
          mm->AddSignalObject(hFit[0]);
          mm->AddArrayObject(n,(TObject**)&hFit[1]);
          mm->AddBackgroundObject(hTotBG);
          mm->AddTotalObject(hTot);
          mm->AddFitObject(fGauss);
          mm->Write();

          // add mean and sigma cut positions
          gMean->SetPoint(j,hData->GetXaxis()->GetBinCenter(j+1),fGauss->GetParameter(1));
          gSigma->SetPoint(j,hData->GetXaxis()->GetBinCenter(j+1),fGauss->GetParameter(2));
          gMean->SetPointError(j,0.,0.);
          gSigma->SetPointError(j,0.,0.);

          delete [] hFit;
          delete fFit;
          delete fGauss;
          delete hDataP;
          delete hTot;
       }

       delete [] hBG;
       delete hData;
       delete hMC;

       // write mean and sigma cut positions
       gMean->Write(Form("Mean_%d", i));
       cut->SetMeanGraph(i, gMean);
       gSigma->Write(Form("Sigma_%d", i));
       cut->SetSigmaGraph(i, gSigma);

       delete gMean;
       delete gSigma;

    }
    delete [] fBG;

    // write cut positions
    cut->Write();
    delete cut;
}
//______________________________________________________________________________
void TMCuts::FitCop(TMObjectCollection* m, const Char_t* s, Int_t iCT, Int_t iRebinMM, Bool_t kFit, Double_t fLeft, Double_t fRight)
{
    Char_t szName[256];
    Char_t szzName[256];

    TFile* fData = 0;
    TFile* fMC = 0;
    TFile* fMC1 = 0;
    TFile* fMC2 = 0;
    TFile** fBG = 0;
    Int_t n = m->GetArraySize();
    fBG = new TFile*[n];

    TIter next(m->GetCollection());
    TObject* o;
    while ((o = (TObject*)next()))
    {
       if (TMTools::IsObjectType(o,"TFile*"))
       {
          Error("TMCuts::FitCop","Collection Objects must be of type TFile");
          gSystem->Exit(0);
       }
    }

    if (!(fData=(TFile*)m->GetDataObject()))
       Error("TMCuts::FitCop","No Data Object found");
    else
       Info("TMCuts::FitCop","Found Data Object %s",fData->GetName());
    if (!(fMC=(TFile*)m->GetSignalObject()))
    {
       if (!(fMC1=(TFile*)m->GetSignal1Object()))
          Error("TMCuts::FitCop","No MC Signal1 Object found");
       else
          Info("TMCuts::FitCop","Found MC Signal1 Object %s",fMC1->GetName());
       if (!(fMC2=(TFile*)m->GetSignal2Object()))
          Error("TMCuts::FitCop","No MC Signal2 Object found");
       else
          Info("TMCuts::FitCop","Found MC Signal2 Object %s",fMC2->GetName());
    }
    else
       Info("TMCuts::FitCop","Found MC Signal Object %s",fMC->GetName());
    if (!(fBG[0]=(TFile*)m->GetBackgroundObject()) && !(fBG=(TFile**)m->GetArrayObject()))
       Error("TMCuts::FitCop","No MC Background Objects found");
    else
       for (Int_t i = 0; i < n; i++)
          Info("TMCuts::FitCop","Found MC Background Object %s",fBG[i]->GetName());

    Info("TMCuts::FitCop","Using %i CT Bins",iCT);

    TOEnergyThetaCut* cut = 0;

    if (kFit) cut = new TOEnergyThetaCut(iCT, "cop_cut");

    for (Int_t i = 0; i < iCT; i++)
    {
       sprintf(szzName,"%s_%i",s,i);
       if (!TMTools::CheckObject(szzName,fData))
          Error("TMCuts::FitCop","Histogram %s in file %s not found", szzName,fData->GetName());
       else
          Info("TMCuts::FitCop","Found Histogram %s in file %s", szzName,fData->GetName());
       if (fMC)
       {
          if (!TMTools::CheckObject(szzName,fMC))
             Error("TMCuts::FitCop","Histogram %s in file %s not found", szzName,fMC->GetName());
          else
             Info("TMCuts::FitCop","Found Histogram %s in file %s", szzName,fMC->GetName());
       }
       else if (fMC1 && fMC2)
       {
          if (!TMTools::CheckObject(szzName,fMC1))
             Error("TMCuts::FitCop","Histogram %s in file %s not found", szzName,fMC1->GetName());
          else
             Info("TMCuts::FitCop","Found Histogram %s in file %s", szzName,fMC1->GetName());
          if (!TMTools::CheckObject(szzName,fMC2))
             Error("TMCuts::FitCop","Histogram %s in file %s not found", szzName,fMC2->GetName());
          else
             Info("TMCuts::FitCop","Found Histogram %s in file %s", szzName,fMC2->GetName());
       }
       else
          Error("TMCuts::FitCop","No Signal Histograms found");
       for (Int_t j = 0; j < n; j++)
       {
          if (!TMTools::CheckObject(szzName,fBG[j]))
             Error("TMCuts::FitCop","Histogram %s in file %s not found", szzName,fBG[j]->GetName());
          else
             Info("TMCuts::FitCop","Found Histogram %s in file %s", szzName,fBG[j]->GetName());
       }

       TH2F* hData = (TH2F*)fData->Get(szzName)->Clone("hData"); 
       TH2F* hMC = 0;
       TH2F* hMC1 = 0;
       TH2F* hMC2 = 0;

       if (fMC)
          hMC = (TH2F*)fMC->Get(szzName)->Clone("hMC"); 
       else 
       {
          hMC1 = (TH2F*)fMC1->Get(szzName)->Clone("hMC1"); 
          hMC2 = (TH2F*)fMC2->Get(szzName)->Clone("hMC2"); 
       }
       TH2F** hBG;
       hBG = new TH2F*[n];
    
       for (Int_t j = 0; j < n; j++)
       {
          sprintf(szName,"hBG_%i",j);
          hBG[j] = (TH2F*)fBG[j]->Get(szzName)->Clone(szName);
       }

       TGraphErrors* gMean = 0;
       TGraphErrors* gSigma = 0;

       if (kFit)
       {
          gMean  = new TGraphErrors(hData->GetNbinsX());
          gSigma = new TGraphErrors(hData->GetNbinsX());
       }

       for (Int_t j = 0; j < hData->GetNbinsX(); j++)
       {
          sprintf(szName,"hDataP_%i_%i",i+1,j+1);
          TH1D* hDataP = (TH1D*)hData->ProjectionY(szName,j+1,j+1,"e");
          hDataP->Rebin(iRebinMM);

          TH1** hFit;
          hFit = new TH1*[n+1];

          TF1* fFit = 0;
          Int_t nMax = 0;

          if (hMC)
          {
             sprintf(szName,"hMCP_%i_%i",i+1,j+1);
             hFit[0] = (TH1D*)hMC->ProjectionY(szName,j+1,j+1,"e");
             hFit[0]->Rebin(iRebinMM);
             for (Int_t k = 0; k < n; k++)
             {
                sprintf(szName,"hBGP_%i_%i_%i",k,i+1,j+1);
                hFit[k+1] = (TH1D*)hBG[k]->ProjectionY(szName,j+1,j+1,"e");
                hFit[k+1]->Rebin(iRebinMM);
             }
   
             TMFunctor::InitFunctor(n+1,hFit);
             fFit = new TF1("fFit",TMFunctor::MyFitFuncMC,hDataP->GetXaxis()->GetBinLowEdge(1),hDataP->GetXaxis()->GetBinUpEdge(hDataP->GetNbinsX()),n+1);
             Double_t pFit[n+1];
             for (Int_t k = 0; k < n+1; k++)
                pFit[k]=1.;
             fFit->SetParameters(pFit);
             for (Int_t k = 0; k < n+1; k++)
             {
                fFit->SetParLimits(k,0.,100000.);
             }
             nMax = n+1;
          }
          else
          {
             sprintf(szName,"hMC1P_%i_%i",i+1,j+1);
             hFit[0] = (TH1D*)hMC1->ProjectionY(szName,j+1,j+1,"e");
             hFit[0]->Rebin(iRebinMM);
             sprintf(szName,"hMC2P_%i_%i",i+1,j+1);
             hFit[1] = (TH1D*)hMC2->ProjectionY(szName,j+1,j+1,"e");
             hFit[1]->Rebin(iRebinMM);

             for (Int_t k = 0; k < n; k++)
             {
                sprintf(szName,"hBGP_%i_%i_%i",k,i+1,j+1);
                hFit[k+2] = (TH1D*)hBG[k]->ProjectionY(szName,j+1,j+1,"e");
                hFit[k+2]->Rebin(iRebinMM);
             }
   
             TMFunctor::InitFunctor(n+2,hFit);
             fFit = new TF1("fFit",TMFunctor::MyFitFuncMC,hDataP->GetXaxis()->GetBinLowEdge(1),hDataP->GetXaxis()->GetBinUpEdge(hDataP->GetNbinsX()),n+2);
             Double_t pFit[n+2];
		pFit[0] = hDataP->GetMaximum();
             for (Int_t k = 1; k < n+2; k++)
                pFit[k]=pFit[0]*0.01;
//             Double_t pFit[n+2];
//             for (Int_t k = 0; k < n+2; k++)
//                pFit[k]=10000.;
             fFit->SetParameters(pFit);
             for (Int_t k = 0; k < n+2; k++)
             {
                fFit->SetParLimits(k,0.,100000.);
             }
             nMax = n+2;
          }

          hDataP->Fit(fFit,"+RW0MQ");

          hFit[0]->Scale(fFit->GetParameter(0));
          sprintf(szName,"hTot_%i_%i",i+1,j+1);
          TH1D* hTot = (TH1D*)hFit[0]->Clone(szName);
//printf("Height Bin 80: %s: %5.3f  factor: %5.3f\n", hFit[0]->GetName(),hFit[0]->GetBinContent(72),fFit->GetParameter(0));
//printf("Height Bin 80: %s: %5.3f\n", hTot->GetName(),hTot->GetBinContent(72));
          for (Int_t k = 1; k < nMax; k++)
          {
             hFit[k]->Scale(fFit->GetParameter(k));
//printf("Height Bin 80: %s: %5.3f  factor: %5.3f\n", hFit[k]->GetName(),hFit[k]->GetBinContent(72),fFit->GetParameter(k));
             hTot->Add(hFit[k]);
//printf("Height Bin 80: %s: %5.3f\n", hTot->GetName(),hTot->GetBinContent(72));

          }

          sprintf(szName,"copObj_%i_%i",i+1,j+1);
          sprintf(szzName,"E_{#gamma}=%5.3f MeV",hData->GetXaxis()->GetBinCenter(j+1));
          TMObjectCollection* mm = new TMObjectCollection(szName,szzName);
          mm->AddDataObject(hDataP);

          TH1* hSig = 0;
          if (hMC)
          { 
             mm->AddSignalObject(hFit[0]);
             mm->AddArrayObject(n,(TObject**)&hFit[1]);
          }
          else
          {
             sprintf(szName,"hSig_%i_%i",i+1,j+1);
             hSig = (TH1*)hFit[0]->Clone(szName);
             hSig->Add(hFit[1]);
             mm->AddSignalObject(hSig);
             mm->AddArrayObject(n,(TObject**)&hFit[2]);
          }
          mm->AddTotalObject(hTot);

          TF1* fGauss = 0;
          if (kFit)
          {
             sprintf(szName,"fGauss_%i_%i",i+1,j+1);
             fGauss = new TF1(szName,
                                   "gaus",
                                   130.,
                                   230.);
             sprintf(szName,"%5.3f",hDataP->GetXaxis()->GetBinCenter(j+1));
             fGauss->SetTitle(szName);
             fGauss->SetParameters(1000.,0.,50.);
             fGauss->SetParLimits(0,0.7*hFit[0]->GetMaximum(),1.3*hFit[0]->GetMaximum());
             fGauss->SetParLimits(1,170.,190.);
             fGauss->SetParLimits(2,1.,100000.);
//             hTot->Fit(fGauss,"+RB0MQ","",fLeft,fRight);
             hFit[0]->Fit(fGauss,"+RB0MQ","",fLeft,fRight);
//             hFit[1]->Fit(fGauss,"+RB0MQ","",fLeft,fRight);
//             if (hMC)
//                hFit[0]->Fit(fGauss,"+RB0MQ");
//             else
//                hSig->Fit(fGauss,"+RB0MQ");
             mm->AddFitObject(fGauss);
             gMean->SetPoint(j,hData->GetXaxis()->GetBinCenter(j+1),fGauss->GetParameter(1));
             gSigma->SetPoint(j,hData->GetXaxis()->GetBinCenter(j+1),fGauss->GetParameter(2));
             gMean->SetPointError(j,0.,0.);
             gSigma->SetPointError(j,0.,0.);
//             gMean->SetPointError(j,hData->GetXaxis()->GetBinWidth(i+1)/2.,fGauss->GetParError(1));
//             gSigma->SetPointError(j,hData->GetXaxis()->GetBinWidth(i+1)/2.,fGauss->GetParError(2));
          }


          mm->Write();

          delete hDataP;
          delete hTot;
          delete [] hFit;
          delete fFit;
          if (kFit) delete fGauss;
       }
       delete [] hBG;
       delete hData;
       if (hMC)
          delete hMC;
       else
       {
          delete hMC1;
          delete hMC2;
       }

       if (kFit)
       {  
//          TOHUtils::RebinGraph(gMean, 2);
//          TOHUtils::SmoothGraph(gMean);
          sprintf(szName, "Mean_%d", i);
          gMean->Write(szName);
          cut->SetMeanGraph(i, gMean);
//          TOHUtils::RebinGraph(gSigma, 2);
//          TOHUtils::SmoothGraph(gSigma);
          sprintf(szName, "Sigma_%d", i);
          gSigma->Write(szName);
          cut->SetSigmaGraph(i, gSigma);

          delete gMean;
          delete gSigma;
       }
    }
    delete [] fBG;

    if (kFit)
    {
       cut->Write();
       delete cut;
    }
}
//______________________________________________________________________________
void TMCuts::FitDEvsE(TH2F* h, Int_t iR, Double_t fSigma, Bool_t kMC)
{
   Char_t szName[256];
   Char_t szFile[256] = "dE_E.cfg";

   ifstream fin;
   TString strLine;
   fin.open(szFile);
   if( !fin.is_open() )
   {
      Error("TMCuts::FitDEvsE","Could not open file %s!", szFile);
      gSystem->Exit(0);
   }

   Int_t nBinsX = h->GetNbinsX();
   Int_t nSlices = nBinsX/iR;
   if (nBinsX%iR)
   {
      nSlices += 1;
      Warning("TMCuts::FitDEvsE","Number of bins %i in histogram %s not a multiple of %i", 
              nBinsX, h->GetName(), iR);
   }

   TGraphErrors* gProtonM = new TGraphErrors(nSlices);
   gProtonM->SetMarkerColor(2);
   gProtonM->SetMarkerStyle(20);
   TGraphErrors* gProtonSP = new TGraphErrors(nSlices);
   gProtonSP->SetMarkerColor(2);
   gProtonSP->SetMarkerStyle(4);
   TGraphErrors* gProtonSM = new TGraphErrors(nSlices);
   gProtonSM->SetMarkerColor(2);
   gProtonSM->SetMarkerStyle(4);

   TGraphErrors* gPionM = new TGraphErrors(nSlices);
   gPionM->SetMarkerColor(6);
   gPionM->SetMarkerStyle(20);
   TGraphErrors* gPionSP = new TGraphErrors(nSlices);
   gPionSP->SetMarkerColor(6);
   gPionSP->SetMarkerStyle(4);
   TGraphErrors* gPionSM = new TGraphErrors(nSlices);
   gPionSM->SetMarkerColor(6);
   gPionSM->SetMarkerStyle(4);

   TGraphErrors* gDeuteronM = 0;
   TGraphErrors* gDeuteronSP = 0;
   TGraphErrors* gDeuteronSM = 0;

   if (!kMC)
   {
      gDeuteronM = new TGraphErrors(nSlices);
      gDeuteronM->SetMarkerColor(4);
      gDeuteronM->SetMarkerStyle(20);
      gDeuteronSP = new TGraphErrors(nSlices);
      gDeuteronSP->SetMarkerColor(4);
      gDeuteronSP->SetMarkerStyle(4);
      gDeuteronSM = new TGraphErrors(nSlices);
      gDeuteronSM->SetMarkerColor(4);
      gDeuteronSM->SetMarkerStyle(4);
   }

   for (Int_t i = 0; i < nSlices; i++) 
   {
      TH1D* hp = 0;
      Double_t pos = 0.;
      if (i < nSlices-1)
      {
         hp = (TH1D*)h->ProjectionY("_py",i*iR+1,i*iR+iR,"e");
         pos = h->GetXaxis()->GetBinLowEdge(i*iR+1)+
                  0.5*(h->GetXaxis()->GetBinUpEdge(i*iR+iR)-
                       h->GetXaxis()->GetBinLowEdge(i*iR+1));
      }
      else
      {
         hp = (TH1D*)h->ProjectionY("_py",i*iR+1,nBinsX,"e");
         pos = h->GetXaxis()->GetBinLowEdge(i*iR+1)+
                  0.5*(h->GetXaxis()->GetBinUpEdge(nBinsX)-
                       h->GetXaxis()->GetBinLowEdge(i*iR+1));
      }

      strLine.ReadLine(fin);
      Double_t ma = 0.;
      Double_t mb = 0.;
      Double_t sa = 0.;
      Double_t sb = 0.;
      Double_t mc;
      Double_t sc;
      if (!kMC)
      {
         mc = 0.;
         sc = 0.;
      }
      Double_t lr = 0.;
      Double_t rr = 0.;

      if ( !kMC && !strLine.BeginsWith("#") && !strLine.BeginsWith(" ") )
         sscanf( strLine.Data(), "%*i %lf %lf %lf %lf %lf %lf %lf %lf",&lr,&rr,&ma,&sa,&mb,&sb,&mc,&sc);
      else
         sscanf( strLine.Data(), "%*i %lf %lf %lf %lf %lf %lf",&lr,&rr,&ma,&sa,&mb,&sb);

      TF1* f = 0;
      if (!kMC)
      {
         f = new TF1("f","pol2(0)+landau(3)+gaus(6)+gaus(9)");
         f->SetLineColor(2);
         f->SetParameters(0,1.);
         f->SetParameters(1,1.);
         f->SetParameters(2,1.);
         f->SetParameters(3,1000.);
         f->SetParameters(4,ma);
         f->SetParameters(5,2.);
         f->SetParameters(6,1000.);
         f->SetParameters(7,mb);
         f->SetParameters(8,2.);
         f->SetParameters(9,1000.);
         f->SetParameters(10,mc);
         f->SetParameters(11,2.);
   
   //      f->SetParLimits(3,50.,10000.);
         f->SetParLimits(4,0.8*ma,1.2*ma);
         f->SetParLimits(5,0.1,5.);
   //      f->SetParLimits(6,50.,10000.);
         f->SetParLimits(7,0.9*mb,1.1*mb);
         f->SetParLimits(8,0.1,5.);
   //      f->SetParLimits(9,50.,10000.);
         f->SetParLimits(10,0.9*mc,1.1*mc);
         f->SetParLimits(11,0.1,5.);
      }
      else
      {
         f = new TF1("f","pol2(0)+gaus(3)");
         f->SetParameter(0,1.);
         f->SetParameter(1,1.);
         f->SetParameter(2,1.);
//         f->SetParameters(3,1.);
//         f->SetParameter(3,1000.);
//         f->SetParameter(4,mb);
//         f->SetParameter(5,sb);

         f->SetParameter(3,hp->GetMaximum());
         f->SetParameter(4,hp->GetBinCenter(hp->GetMaximumBin()));
         f->SetParameter(5,0.5);

         f->SetParLimits(3,0.8*hp->GetMaximum(),1.2*hp->GetMaximum());
         f->SetParLimits(4,0.9*hp->GetBinCenter(hp->GetMaximumBin()),1.1*hp->GetBinCenter(hp->GetMaximumBin()));
         f->SetParLimits(5,0.9*0.5,1.1*0.5);

   //      f->SetParLimits(3,50.,10000.);
//         f->SetParLimits(4,0.9*mb,1.1*mb);
//         f->SetParLimits(5,0.9*sb,1.1*sb);
      }
      f->SetLineColor(2);
      f->SetRange(lr,rr);

      hp->Fit(f,"+RB0Q");
      if (!kMC)
         Info("TMCuts::FitDEvsE","%i - range:%5.3f-%5.3f - mean1: %5.3f sigma1: %5.3f   mean2: %5.3f sigma2: %5.3f   mean3: %5.3f sigma3: %5.3f\n", i, lr,rr,f->GetParameter(4),f->GetParameter(5),f->GetParameter(7),f->GetParameter(8),f->GetParameter(10),f->GetParameter(11));
      else
         Info("TMCuts::FitDEvsE","%i - range:%5.3f-%5.3f - mean1: %5.3f sigma1: %5.3f\n", i, lr,rr,f->GetParameter(4),f->GetParameter(5));

      if (!kMC)
      {
         gPionM->SetPoint(i,pos,f->GetParameter(4));
         gProtonM->SetPoint(i,pos,f->GetParameter(7));
         gPionM->SetPointError(i,0.,f->GetParError(4));
         gProtonM->SetPointError(i,0.,f->GetParError(7));
         gPionSM->SetPoint(i,pos,f->GetParameter(4)-fSigma*f->GetParameter(5));
         gProtonSM->SetPoint(i,pos,f->GetParameter(7)-fSigma*f->GetParameter(8));
         gPionSM->SetPointError(i,0.,f->GetParError(4)+fSigma*f->GetParError(5));
         gProtonSM->SetPointError(i,0.,f->GetParError(7)+fSigma*f->GetParError(8));
         gPionSP->SetPoint(i,pos,f->GetParameter(4)+fSigma*f->GetParameter(5));
         gProtonSP->SetPoint(i,pos,f->GetParameter(7)+fSigma*f->GetParameter(8));
         gPionSP->SetPointError(i,0.,f->GetParError(4)+fSigma*f->GetParError(5));
         gProtonSP->SetPointError(i,0.,f->GetParError(7)+fSigma*f->GetParError(8));

//         gDeuteronM->SetPoint(i,pos,f->GetParameter(10));
//         gDeuteronM->SetPointError(i,0.,f->GetParError(10));
//         gDeuteronSM->SetPoint(i,pos,f->GetParameter(10)-f->GetParameter(11));
//         gDeuteronSM->SetPointError(i,0.,f->GetParError(10)+f->GetParError(11));
//         gDeuteronSP->SetPoint(i,pos,f->GetParameter(10)+f->GetParameter(11));
//         gDeuteronSP->SetPointError(i,0.,f->GetParError(10)+f->GetParError(11));
      }
      else
      {
         gProtonM->SetPoint(i,pos,f->GetParameter(4));
         gProtonM->SetPointError(i,0.,f->GetParError(4));
         gProtonSM->SetPoint(i,pos,f->GetParameter(4)-fSigma*f->GetParameter(5));
         gProtonSM->SetPointError(i,0.,f->GetParError(4)+fSigma*f->GetParError(5));
         gProtonSP->SetPoint(i,pos,f->GetParameter(4)+fSigma*f->GetParameter(5));
         gProtonSP->SetPointError(i,0.,f->GetParError(4)+fSigma*f->GetParError(5));

         if (lr == 0 && rr == 0)
         {
            gProtonM->RemovePoint(i);
            gProtonSM->RemovePoint(i);
            gProtonSP->RemovePoint(i);
printf("removed\n");
            continue;
         }
      }

      sprintf(szName,"dE_E_%i",i);
      TCanvas* c = new TCanvas(szName,szName,800,600);

      c->cd();
      hp->Draw();
      f->Draw("same");

      c->Write();

      delete hp;
      delete f;
      delete c;
   }

   gProtonM->SetName("proton_mean");
   gProtonSM->SetName("proton_ms");
   gProtonSP->SetName("proton_ps");

   if (!kMC)
   {
      gPionM->SetName("pion_mean");
//      gDeuteronM->SetName("deuteron_mean");
      gPionSM->SetName("pion_ms");
//      gDeuteronSM->SetName("deuteron_m1s");
      gPionSP->SetName("pion_ps");
//      gDeuteronSP->SetName("deuteron_p1s");
   }

//   TF1* fPionM      = new TF1("fPionM","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fPionM->SetParameters(1.,1.,0.5);
//   TF1* fProtonM    = new TF1("fProtonM","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fProtonM->SetParameters(1.,1.,0.5);
//   TF1* fDeuteronM  = new TF1("fDeuteronM","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fDeuteronM->SetParameters(1.,1.,0.5);
//   TF1* fPionSM     = new TF1("fPionSM","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fPionSM->SetParameters(1.,1.,0.5);
//   TF1* fProtonSM   = new TF1("fProtonSM","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fProtonSM->SetParameters(1.,1.,0.5);
//   TF1* fDeuteronSM = new TF1("fDeuteronSM","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fDeuteronSM->SetParameters(1.,1.,0.5);
//   TF1* fPionSP     = new TF1("fPionSP","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fPionSP->SetParameters(1.,1.,0.5);
//   TF1* fProtonSP   = new TF1("fProtonSP","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fProtonSP->SetParameters(1.,1.,0.5);
//   TF1* fDeuteronSP = new TF1("fDeuteronSP","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fDeuteronSP->SetParameters(1.,1.,0.5);


   TF1* fProtonM    = new TF1("fProtonM",   "pol5",0.,1000.);
   TF1* fProtonSM   = new TF1("fProtonSM",  "pol5",0.,1000.);
   TF1* fProtonSP   = new TF1("fProtonSP",  "pol5",0.,1000.);

   Double_t par[6];
   fProtonM->SetParameters(1.,1.,1.,1.,1.,1.);
   fProtonSM->SetParameters(1.,1.,1.,1.,1.,1.);
   fProtonSP->SetParameters(1.,1.,1.,1.,1.,1.);

   fProtonM->SetLineColor(2);
   fProtonSM->SetLineColor(2);
   fProtonSP->SetLineColor(2);
   fProtonSM->SetLineStyle(2);
   fProtonSP->SetLineStyle(2);

   TF1* fPionM = 0;
   TF1* fPionSM = 0;
   TF1* fPionSP = 0;
//   TF1* fDeuteronM = 0;
//   TF1* fDeuteronSM = 0;
//   TF1* fDeuteronSP = 0;

   if (!kMC)
   {
      fPionM      = new TF1("fPionM",    "pol2",0.,1000.);
      fPionSM     = new TF1("fPionSM",    "pol2",0.,1000.);
      fPionSP     = new TF1("fPionSP",    "pol2",0.,1000.);
                                                            
//      fDeuteronM  = new TF1("fDeuteronM", "pol5",0.,1000.);
//      fDeuteronSM = new TF1("fDeuteronSM","pol5",0.,1000.);
//      fDeuteronSP = new TF1("fDeuteronSP","pol5",0.,1000.);

      fPionM->SetLineColor(6);
//      fDeuteronM->SetLineColor(4);
      fPionSM->SetLineColor(6);
//      fDeuteronSM->SetLineColor(4);
      fPionSP->SetLineColor(6);
//      fDeuteronSP->SetLineColor(4);
      fPionSM->SetLineStyle(2);
//      fDeuteronSM->SetLineStyle(2);
      fPionSP->SetLineStyle(2);
//      fDeuteronSP->SetLineStyle(2);
   }

   Double_t fLeft = 0.;
   Double_t fRight = 1000.;

   gProtonM->Fit(fProtonM,"+M0Q","",fLeft,fRight);
   gProtonM->Write();
   fProtonM->SetName("fit_proton_mean");
   fProtonM->Write();

   fProtonM->GetParameters(&par[0]);
   fProtonSM->SetParameters(par);
   fProtonSP->SetParameters(par);

   gProtonSM->Fit(fProtonSM,"+M0Q","",fLeft,fRight);
   gProtonSM->Write();
   fProtonSM->SetName("fit_proton_m1s");
   fProtonSM->Write();

   gProtonSP->Fit(fProtonSP,"+M0Q","",fLeft,fRight);
   gProtonSP->Write();
   fProtonSP->SetName("fit_proton_1s");
   fProtonSP->Write();

   if (!kMC)
   {
      gPionM->Fit(fPionM,"+ME0Q","",fLeft,fRight);
      gPionM->Write();
      fPionM->SetName("fit_pion_mean");
      fPionM->Write();

      fPionM->GetParameters(&par[0]);
      fPionSM->SetParameters(par);
      fPionSP->SetParameters(par);
     
      gPionSM->Fit(fPionSM,"+ME0Q","",fLeft,fRight);
      gPionSM->Write();
      fPionSM->SetName("fit_pion_m1s");
      fPionSM->Write();

      gPionSP->Fit(fPionSP,"+ME0Q","",fLeft,fRight);
      gPionSP->Write();
      fPionSP->SetName("fit_pion_1s");
      fPionSP->Write();
   
//      gDeuteronM->Fit(fDeuteronM,"+ME0Q","",fLeft,fRight);
//      gDeuteronM->Write();
//      fDeuteronM->SetName("fit_deuteron_mean");
//      fDeuteronM->Write();

//      fDeuteronM->GetParameters(&par[0]);
//      fDeuteronSM->SetParameters(par);
//      fDeuteronSP->SetParameters(par);

//      gDeuteronSM->Fit(fDeuteronSM,"+ME0Q","",fLeft,fRight);
//      gDeuteronSM->Write();
//      fDeuteronSM->SetName("fit_deuteron_m1s");
//      fDeuteronSM->Write();
     
//      gDeuteronSP->Fit(fDeuteronSP,"+ME0Q","",fLeft,fRight);
//      gDeuteronSP->Write();
//      fDeuteronSP->SetName("fit_deuteron_1s");
//      fDeuteronSP->Write();
   }

   TCanvas* d = new TCanvas("d","dE_E",800,600);
   d->cd();

//   gPionM->GetXaxis()->SetRangeUser(0.,1000.);
//   gPionM->GetYaxis()->SetRangeUser(0.,10.);
//   gPionM->Draw("AP");
   gProtonM->GetXaxis()->SetRangeUser(0.,1000.);
   gProtonM->GetYaxis()->SetRangeUser(0.,10.);
   gProtonM->Draw("AP");
//   gDeuteronM->Draw("Psame");
//   gPionSP->Draw("Psame");
   gProtonSP->Draw("Psame");
//   gDeuteronSP->Draw("Psame");
//   gPionSM->Draw("Psame");
   gProtonSM->Draw("Psame");
//   gDeuteronSM->Draw("Psame");

//   fPionM->Draw("same");
   fProtonM->Draw("same");
//   fDeuteronM->Draw("same");
//   fPionSP->Draw("same");
   fProtonSP->Draw("same");
//   fDeuteronSP->Draw("same");
//   fPionSM->Draw("same");
   fProtonSM->Draw("same");
//   fDeuteronSM->Draw("same");

   d->Write();
   delete d;

   TCanvas* e = new TCanvas("e","dE_E 2D",800,600);
   e->cd();

   h->Draw("col");
   fProtonM->SetLineColor(0);
   fProtonSM->SetLineColor(0);
   fProtonSP->SetLineColor(0);
   fProtonM->Draw("same");
   fProtonSM->Draw("same");
   fProtonSP->Draw("same");
   e->Write();

   delete e;

   delete gProtonM;
   delete gProtonSP;
   delete gProtonSM;
   delete fProtonM;
   delete fProtonSP;
   delete fProtonSM;


   delete fPionM;
//   delete fDeuteronM;
   delete fPionSP;
//   delete fDeuteronSP;
   delete fPionSM;
//   delete fDeuteronSM;
   delete gPionM;
//   delete gDeuteronM;
   delete gPionSP;
//   delete gDeuteronSP;
   delete gPionSM;
//   delete gDeuteronSM;

   fin.close();

   return;
}
//______________________________________________________________________________
void TMCuts::FitDEvsTOFx(TH2F* h, Int_t iR, Double_t fSigma, Bool_t kMC)
{
   Char_t szName[256];
   Char_t szFile[256] = "dE_TOF_x.cfg";

   ifstream fin;
   TString strLine;
   fin.open(szFile);
   if( !fin.is_open() )
   {
      Error("TMCuts::FitDEvsTOFx","Could not open file %s!", szFile);
      gSystem->Exit(0);
   }

   Int_t nBinsY = h->GetNbinsY();
   Int_t nSlices = nBinsY/iR;
   if (nBinsY%iR)
   {
      nSlices += 1;
      Warning("TMCuts::FitDEvsTOFx","Number of bins %i in histogram %s not a multiple of %i", 
              nBinsY, h->GetName(), iR);
   }

   TGraphErrors* gProtonM = new TGraphErrors(nSlices);
   gProtonM->SetMarkerColor(2);
   gProtonM->SetMarkerStyle(20);
   TGraphErrors* gProtonSP = new TGraphErrors(nSlices);
   gProtonSP->SetMarkerColor(2);
   gProtonSP->SetMarkerStyle(4);
   TGraphErrors* gProtonSM = new TGraphErrors(nSlices);
   gProtonSM->SetMarkerColor(2);
   gProtonSM->SetMarkerStyle(4);

   TGraphErrors* gPionM = new TGraphErrors(nSlices);
   gPionM->SetMarkerColor(6);
   gPionM->SetMarkerStyle(20);
   TGraphErrors* gPionSP = new TGraphErrors(nSlices);
   gPionSP->SetMarkerColor(6);
   gPionSP->SetMarkerStyle(4);
   TGraphErrors* gPionSM = new TGraphErrors(nSlices);
   gPionSM->SetMarkerColor(6);
   gPionSM->SetMarkerStyle(4);

   TGraphErrors* gDeuteronM = 0;
   TGraphErrors* gDeuteronSP = 0;
   TGraphErrors* gDeuteronSM = 0;

   if (!kMC)
   {
      gDeuteronM = new TGraphErrors(nSlices);
      gDeuteronM->SetMarkerColor(4);
      gDeuteronM->SetMarkerStyle(20);
      gDeuteronSP = new TGraphErrors(nSlices);
      gDeuteronSP->SetMarkerColor(4);
      gDeuteronSP->SetMarkerStyle(4);
      gDeuteronSM = new TGraphErrors(nSlices);
      gDeuteronSM->SetMarkerColor(4);
      gDeuteronSM->SetMarkerStyle(4);
   }

   for (Int_t i = 0; i < nSlices; i++) 
   {
      TH1D* hp = 0;
      Double_t pos = 0.;
      if (i < nSlices-1)
      {
         hp = (TH1D*)h->ProjectionX("_px",i*iR+1,i*iR+iR);
         pos = h->GetYaxis()->GetBinLowEdge(i*iR+1)+
                  0.5*(h->GetYaxis()->GetBinUpEdge(i*iR+iR)-
                       h->GetYaxis()->GetBinLowEdge(i*iR+1));
      }
      else
      {
         hp = (TH1D*)h->ProjectionX("_px",i*iR+1,nBinsY);
         pos = h->GetYaxis()->GetBinLowEdge(i*iR+1)+
                  0.5*(h->GetYaxis()->GetBinUpEdge(nBinsY)-
                       h->GetYaxis()->GetBinLowEdge(i*iR+1));
      }

      strLine.ReadLine(fin);
      Double_t ma = 0.;
      Double_t mb = 0.;
      Double_t sa = 0.;
      Double_t sb = 0.;
      Double_t mc;
      Double_t sc;
      if (!kMC)
      {
         mc = 0.;
         sc = 0.;
      }
      Double_t lr = 0.;
      Double_t rr = 0.;

      if ( !kMC && !strLine.BeginsWith("#") && !strLine.BeginsWith(" ") )
         sscanf( strLine.Data(), "%*i %lf %lf %lf %lf %lf %lf %lf %lf",&lr,&rr,&ma,&sa,&mb,&sb,&mc,&sc);
      else
         sscanf( strLine.Data(), "%*i %lf %lf %lf %lf %lf %lf",&lr,&rr,&ma,&sa,&mb,&sb);

      TF1* f = 0;
      if (!kMC)
      {
         f = new TF1("f","pol2(0)+gaus(3)+gaus(6)");
         f->SetParameters(0,1.);
         f->SetParameters(1,1.);
         f->SetParameters(2,1.);
         f->SetParameters(3,1000.);
         f->SetParameters(4,ma);
         f->SetParameters(5,sa);
         f->SetParameters(6,1000.);
         f->SetParameters(7,mb);
         f->SetParameters(8,sb);
         if (i >= 35) f->FixParameter(3.0,0);
         else f->SetParLimits(3,0.0,10000.);
         f->SetParLimits(4,0.8*ma,1.2*ma);
         f->SetParLimits(5,0.9*sa,1.1*sa);
         f->SetParLimits(6,0.0,10000.);
         f->SetParLimits(7,0.9*mb,1.1*mb);
         f->SetParLimits(8,0.8*sa,1.2*sb);

//         f->SetParameters(9,1000.);
//         f->SetParameters(10,mc);
//         f->SetParameters(11,2.);

//         f->SetParLimits(9,50.,10000.);
//         f->SetParLimits(10,0.95*mc,1.15*mc);
//         f->SetParLimits(11,0.9*sa,1.1*sc);
      }
      else
      {
         f = new TF1("f","pol2(0)+gaus(3)");
         f->SetParameters(0,1.);
         f->SetParameters(1,1.);
         f->SetParameters(2,1.);
//         f->SetParameters(3,1.);
         f->SetParameters(3,1000.);
         f->SetParameters(4,mb);
         f->SetParameters(5,sb);

   //      f->SetParLimits(3,50.,10000.);
         f->SetParLimits(4,0.9*mb,1.1*mb);
         f->SetParLimits(5,0.9*sb,1.1*sb);
      }
      f->SetLineColor(2);
//      f->SetRange(h->GetXaxis()->GetBinLowEdge(1),h->GetXaxis()->GetBinUpEdge(nBinsY));
      f->SetRange(lr,rr);

      hp->Fit(f,"+RMB0Q");
      if (!kMC)
         Info("TMCuts::FitDEvsTOFx","%i - range:%5.3f-%5.3f - mean1: %5.3f sigma1: %5.3f   mean2: %5.3f sigma2: %5.3f\n", i, lr,rr,f->GetParameter(4),f->GetParameter(5),f->GetParameter(7),f->GetParameter(8));
      else
         Info("TMCuts::FitDEvsTOFx","%i - range:%5.3f-%5.3f - mean1: %5.3f sigma1: %5.3f\n", i, lr,rr,f->GetParameter(4),f->GetParameter(5));

      if (!kMC)
      {
         gProtonM->SetPoint(i,f->GetParameter(7),pos);
         gProtonM->SetPointError(i,f->GetParError(7),0.);
         gProtonSM->SetPoint(i,f->GetParameter(7)-fSigma*f->GetParameter(8),pos);
         gProtonSM->SetPointError(i,f->GetParError(7)+fSigma*f->GetParError(8),0.);
         gProtonSP->SetPoint(i,f->GetParameter(7)+fSigma*f->GetParameter(8),pos);
         gProtonSP->SetPointError(i,f->GetParError(7)+fSigma*f->GetParError(8),0.);
   
         gPionM->SetPoint(i,f->GetParameter(4),pos);
         gPionM->SetPointError(i,f->GetParError(4),0.);
         gPionSM->SetPoint(i,f->GetParameter(4)-f->GetParameter(5),pos);
         gPionSM->SetPointError(i,f->GetParError(4)+f->GetParError(5),0.);
         gPionSP->SetPoint(i,f->GetParameter(4)+f->GetParameter(5),pos);
         gPionSP->SetPointError(i,f->GetParError(4)+f->GetParError(5),0.);

//         gDeuteronM->SetPoint(i,f->GetParameter(10),pos);
//         gDeuteronM->SetPointError(i,f->GetParError(10),0.);
//         gDeuteronSM->SetPoint(i,f->GetParameter(10)-f->GetParameter(11),pos);
//         gDeuteronSM->SetPointError(i,f->GetParError(10)+f->GetParError(11),0.);
//         gDeuteronSP->SetPoint(i,f->GetParameter(10)+f->GetParameter(11),pos);
//         gDeuteronSP->SetPointError(i,f->GetParError(10)+f->GetParError(11),0.);
      }
      else
      {
         gProtonM->SetPoint(i,f->GetParameter(4),pos);
         gProtonM->SetPointError(i,f->GetParError(4),0.);
         gProtonSM->SetPoint(i,f->GetParameter(4)-fSigma*f->GetParameter(5),pos);
         gProtonSM->SetPointError(i,f->GetParError(4)+fSigma*f->GetParError(5),0.);
         gProtonSP->SetPoint(i,f->GetParameter(4)+fSigma*f->GetParameter(5),pos);
         gProtonSP->SetPointError(i,f->GetParError(4)+fSigma*f->GetParError(5),0.);

         if (lr == 0 && rr == 0)
         {
            gProtonM->RemovePoint(i);
            gProtonSM->RemovePoint(i);
            gProtonSP->RemovePoint(i);
printf("removed\n");
            continue;
         }   
      }

      sprintf(szName,"dE_TOF_x_%i",i);
      TCanvas* c = new TCanvas(szName,szName,800,600);

      c->cd();
      hp->Draw();
      hp->GetXaxis()->SetRangeUser(0.,20.);
      f->Draw("same");

      c->Write();

      delete hp;
      delete f;
      delete c;
   }

   gProtonM->SetName("proton_mean_x");
   gProtonSM->SetName("proton_ms_x");
   gProtonSP->SetName("proton_ps_x");
   if (!kMC)
   {
      gPionM->SetName("pion_mean_x");
      gPionSM->SetName("pion_m3s_x");
      gPionSP->SetName("pion_p3s_x");
//      gDeuteronM->SetName("deuteron_mean_x");
//      gDeuteronSM->SetName("deuteron_m3s_x");
//      gDeuteronSP->SetName("deuteron_p3s_x");
   }

//   TF1* fPionM      = new TF1("fPionM","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fPionM->SetParameters(1.,1.,0.5);
//   TF1* fProtonM    = new TF1("fProtonM","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fProtonM->SetParameters(1.,1.,0.5);
//   TF1* fDeuteronM  = new TF1("fDeuteronM","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fDeuteronM->SetParameters(1.,1.,0.5);
//   TF1* fPionSM     = new TF1("fPionSM","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fPionSM->SetParameters(1.,1.,0.5);
//   TF1* fProtonSM   = new TF1("fProtonSM","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fProtonSM->SetParameters(1.,1.,0.5);
//   TF1* fDeuteronSM = new TF1("fDeuteronSM","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fDeuteronSM->SetParameters(1.,1.,0.5);
//   TF1* fPionSP     = new TF1("fPionSP","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fPionSP->SetParameters(1.,1.,0.5);
//   TF1* fProtonSP   = new TF1("fProtonSP","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fProtonSP->SetParameters(1.,1.,0.5);
//   TF1* fDeuteronSP = new TF1("fDeuteronSP","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fDeuteronSP->SetParameters(1.,1.,0.5);

   TF1* fProtonM    = new TF1("fProtonM",   "pol5",0.,20.);//[0]+[1]/TMath::Power(x-[2],[3])",4.,12.);
   TF1* fProtonSM   = new TF1("fProtonSM",  "pol5",0.,20.);
   TF1* fProtonSP   = new TF1("fProtonSP",  "pol5",0.,20.);
   Double_t par[6];
   fProtonM->SetParameters(1.,1.,1.,1.,1.,1.);
   fProtonSM->SetParameters(1.,1.,1.,1.,1.,1.);
   fProtonSP->SetParameters(1.,1.,1.,1.,1.,1.);

   fProtonM->SetLineColor(2);
   fProtonSM->SetLineColor(2);
   fProtonSM->SetLineStyle(2);
   fProtonSP->SetLineColor(2);
   fProtonSP->SetLineStyle(2);

   TF1* fPionM = 0;
   TF1* fPionSM = 0;
   TF1* fPionSP = 0;
//   TF1* fDeuteronM = 0;
//   TF1* fDeuteronSM = 0;
//   TF1* fDeuteronSP = 0;

   if (!kMC)
   {
      fPionM      = new TF1("fPionM",     "pol5",0.,20.);
      fPionSM     = new TF1("fPionSM",    "pol5",0.,20.);
      fPionSP     = new TF1("fPionSP",    "pol5",0.,20.);
   
      fPionM->SetLineColor(6);
      fPionSM->SetLineColor(6);
      fPionSP->SetLineColor(6);
      fPionSM->SetLineStyle(2);
      fPionSP->SetLineStyle(2);

//      fDeuteronM  = new TF1("fDeuteronM", "pol5",0.,20.);
//      fDeuteronSM = new TF1("fDeuteronSM","pol5",0.,20.);
//      fDeuteronSP = new TF1("fDeuteronSP","pol5",0.,20.);
//
//      fDeuteronM->SetLineColor(4);
//      fDeuteronSM->SetLineColor(4);
//      fDeuteronSP->SetLineColor(4);
//      fDeuteronSM->SetLineStyle(2);
//      fDeuteronSP->SetLineStyle(2);
   }

   Double_t fLeft = 0.;
   Double_t fRight = 20.;

   gProtonM->Fit(fProtonM,"+M0Q","",4.0,10.5);
   gProtonM->Write();
   fProtonM->SetName("fit_proton_mean_x");
   fProtonM->Write();

   fProtonM->GetParameters(&par[0]);
   fProtonSM->SetParameters(par);
   fProtonSP->SetParameters(par);

   gProtonSM->Fit(fProtonSM,"+M0Q","",2.5,9.0);
   gProtonSM->Write();
   fProtonSM->SetName("fit_proton_m1s_x");
   fProtonSM->Write();

   gProtonSP->Fit(fProtonSP,"+M0Q","",5.5,11.5);
   gProtonSP->Write();
   fProtonSP->SetName("fit_proton_1s_x");
   fProtonSP->Write();

   if (!kMC)
   {
      gPionM->Fit(fPionM,"+M0Q","",fLeft,fRight);
      gPionM->Write();
      fPionM->SetName("fit_pion_mean_x");
      fPionM->Write();
   
      fPionM->GetParameters(&par[0]);
      fPionSM->SetParameters(par);
      fPionSP->SetParameters(par);
   
      gPionSM->Fit(fPionSM,"+M0Q","",fLeft,fRight);
      gPionSM->Write();
      fPionSM->SetName("fit_pion_m1s_x");
      fPionSM->Write();
   
      gPionSP->Fit(fPionSP,"+M0Q","",fLeft,fRight);
      gPionSP->Write();
      fPionSP->SetName("fit_pion_1s_x");
      fPionSP->Write();

//      gDeuteronM->Fit(fDeuteronM,"+M0Q","",fLeft,fRight);
//      gDeuteronM->Write();
//      fDeuteronM->SetName("fit_deuteron_mean_x");
//      fDeuteronM->Write();
//
//      fDeuteronM->GetParameters(&par[0]);
//      fDeuteronSM->SetParameters(par);
//      fDeuteronSP->SetParameters(par);
//
//      gDeuteronSM->Fit(fDeuteronSM,"+M0Q","",fLeft,fRight);
//      gDeuteronSM->Write();
//      fDeuteronSM->SetName("fit_deuteron_m1s_x");
//      fDeuteronSM->Write();
//
//      gDeuteronSP->Fit(fDeuteronSP,"+M0Q","",fLeft,fRight);
//      gDeuteronSP->Write();
//      fDeuteronSP->SetName("fit_deuteron_1s_x");
//      fDeuteronSP->Write();
   }

   TCanvas* d = new TCanvas("dE_TOF_x","dE_TOF_x",800,600);
   d->cd();

   gProtonM->GetXaxis()->SetRangeUser(-50.,50.);
   gProtonM->GetYaxis()->SetRangeUser(0.,500.);
   gProtonM->Draw("AP");
   gProtonSP->Draw("Psame");
   gProtonSM->Draw("Psame");
   fProtonM->Draw("same");
   fProtonSP->Draw("same");
   fProtonSM->Draw("same");

//   gPionM->Draw("AP");
//   gPionSP->Draw("Psame");
//   gPionSM->Draw("Psame");
//   fPionM->Draw("same");
//   fPionSP->Draw("same");
//   fPionSM->Draw("same");

   if (!kMC)
   {
//      gDeuteronM->Draw("Psame");
//      gDeuteronSP->Draw("Psame");
//      gDeuteronSM->Draw("Psame");
//      fDeuteronM->Draw("same");
//      fDeuteronSP->Draw("same");
//      fDeuteronSM->Draw("same");
   }

   d->Write();
   delete d;

   TCanvas* e = new TCanvas("dE_TOF_2D_x","dE_TOF_2D_x",800,600);
   e->cd();

   h->Draw("col");
   fProtonM->SetLineColor(0);
   fProtonSM->SetLineColor(0);
   fProtonSP->SetLineColor(0);
   fProtonM->Draw("same");
   fProtonSM->Draw("same");
   fProtonSP->Draw("same");
   e->Write();

   delete e;

   delete gProtonM;
   delete gProtonSP;
   delete gProtonSM;
   delete fProtonM;
   delete fProtonSP;
   delete fProtonSM;
   if (!kMC)
   {
      delete gPionM;
      delete gPionSP;
      delete gPionSM;
      delete fPionM;
      delete fPionSP;
      delete fPionSM;
//      delete gDeuteronM;
//      delete gDeuteronSP;
//      delete gDeuteronSM;
//      delete fDeuteronM;
//      delete fDeuteronSP;
//      delete fDeuteronSM;
   }

   fin.close();

   return;
}
//______________________________________________________________________________
void TMCuts::FitDEvsTOFy(TH2F* h, Int_t iR)
{
   Char_t szName[256];
   Char_t szFile[256] = "dE_TOF_y.cfg";

   ifstream fin;
   TString strLine;
   fin.open(szFile);
   if( !fin.is_open() )
   {
      Error("TMCuts::FitDEvsTOFy","Could not open file %s!", szFile);
      gSystem->Exit(0);
   }

   Int_t nBinsX = h->GetNbinsX();
   Int_t nSlices = nBinsX/iR;
   if (nBinsX%iR)
   {
      nSlices += 1;
      Warning("TMCuts::FitDEvsTOFy","Number of bins %i in histogram %s not a multiple of %i", 
              nBinsX, h->GetName(), iR);
   }

   TGraphErrors* gProtonM = new TGraphErrors(nSlices);
   TGraphErrors* gDeuteronM = new TGraphErrors(nSlices);
   gProtonM->SetMarkerColor(2);
   gDeuteronM->SetMarkerColor(4);
   gProtonM->SetMarkerStyle(20);
   gDeuteronM->SetMarkerStyle(20);
   TGraphErrors* gProtonSP = new TGraphErrors(nSlices);
   TGraphErrors* gDeuteronSP = new TGraphErrors(nSlices);
   gProtonSP->SetMarkerColor(2);
   gDeuteronSP->SetMarkerColor(4);
   gProtonSP->SetMarkerStyle(4);
   gDeuteronSP->SetMarkerStyle(4);
   TGraphErrors* gProtonSM = new TGraphErrors(nSlices);
   TGraphErrors* gDeuteronSM = new TGraphErrors(nSlices);
   gProtonSM->SetMarkerColor(2);
   gDeuteronSM->SetMarkerColor(4);
   gProtonSM->SetMarkerStyle(4);
   gDeuteronSM->SetMarkerStyle(4);

   for (Int_t i = 0; i < nSlices; i++) 
   {
      TH1D* hp = 0;
      Double_t pos = 0.;
      if (i < nSlices-1)
      {
         hp = (TH1D*)h->ProjectionY("_py",i*iR+1,i*iR+iR,"e");
         pos = h->GetXaxis()->GetBinLowEdge(i*iR+1)+
                  0.5*(h->GetXaxis()->GetBinUpEdge(i*iR+iR)-
                       h->GetXaxis()->GetBinLowEdge(i*iR+1));
      }
      else
      {
         hp = (TH1D*)h->ProjectionY("_py",i*iR+1,nBinsX,"e");
         pos = h->GetXaxis()->GetBinLowEdge(i*iR+1)+
                  0.5*(h->GetXaxis()->GetBinUpEdge(nBinsX)-
                       h->GetXaxis()->GetBinLowEdge(i*iR+1));
      }

      strLine.ReadLine(fin);
      Double_t ma = 0.;
      Double_t mb = 0.;
      Double_t sa = 0.;
      Double_t sb = 0.;
      Double_t lr = 0.;
      Double_t rr = 0.;

      if ( !strLine.BeginsWith("#") && !strLine.BeginsWith(" ") )
         sscanf( strLine.Data(), "%*i %lf %lf %lf %lf %lf %lf",&lr,&rr,&ma,&sa,&mb,&sb);

      TF1* f = 0;
      f = new TF1("f","pol2(0)+gaus(4)+gaus(7)");
      f->SetLineColor(2);
//      f->SetRange(h->GetXaxis()->GetBinLowEdge(1),h->GetXaxis()->GetBinUpEdge(nBinsY));
      f->SetRange(lr,rr);
      f->SetParameters(0,1.);
      f->SetParameters(1,1.);
      f->SetParameters(2,1.);
      f->SetParameters(3,1.);
      f->SetParameters(4,1000.);
      f->SetParameters(5,ma);
      f->SetParameters(6,2.);
      f->SetParameters(7,1000.);
      f->SetParameters(8,mb);
      f->SetParameters(9,2.);
//      f->SetParameters(9,1000.);
//      f->SetParameters(10,mc);
//      f->SetParameters(11,2.);

//      f->SetParLimits(4,50.,10000.);
//      f->SetParLimits(5,0.9*ma,1.1*ma);
      f->SetParLimits(6,0.9*sa,1.1*sa);
//      f->SetParLimits(7,50.,10000.);
//      f->SetParLimits(8,0.9*mb,1.1*mb);
      f->SetParLimits(9,0.9*sa,1.1*sb);
////      f->SetParLimits(9,50.,10000.);
//      f->SetParLimits(10,0.95*mc,1.15*mc);
//      f->SetParLimits(11,0.9*sa,1.1*sc);

      hp->Fit(f,"+RB0Q");
Info("TMCuts::FitDEvsTOFy","%i - range:%5.3f-%5.3f - mean1: %5.3f sigma1: %5.3f   mean2: %5.3f sigma2: %5.3f\n", i, lr,rr,f->GetParameter(5),f->GetParameter(6),f->GetParameter(8),f->GetParameter(9));

      gProtonM->SetPoint(i,f->GetParameter(5),pos);
      gDeuteronM->SetPoint(i,f->GetParameter(8),pos);
      gProtonM->SetPointError(i,0.,f->GetParError(5));
      gDeuteronM->SetPointError(i,0.,f->GetParError(8));
      gProtonSM->SetPoint(i,f->GetParameter(5)-3.*f->GetParameter(6),pos);
      gDeuteronSM->SetPoint(i,f->GetParameter(8)-f->GetParameter(9),pos);
      gProtonSM->SetPointError(i,f->GetParError(5)+f->GetParError(6),0.);
      gDeuteronSM->SetPointError(i,f->GetParError(8)+f->GetParError(9),0.);
      gProtonSP->SetPoint(i,f->GetParameter(5)+3.*f->GetParameter(6),pos);
      gDeuteronSP->SetPoint(i,f->GetParameter(8)+f->GetParameter(9),pos);
      gProtonSP->SetPointError(i,f->GetParError(5)+f->GetParError(6),0.);
      gDeuteronSP->SetPointError(i,f->GetParError(8)+f->GetParError(9),0.);

      sprintf(szName,"dE_TOF_y_%i",i);
      TCanvas* c = new TCanvas(szName,szName,800,600);

      c->cd();
      hp->Draw();
      hp->GetXaxis()->SetRangeUser(0.,500.);
      f->Draw("same");

      c->Write();

      delete hp;
      delete f;
      delete c;
   }

   gProtonM->SetName("proton_mean_y");
   gDeuteronM->SetName("deuteron_mean_y");
   gProtonSM->SetName("proton_m1s_y");
   gDeuteronSM->SetName("deuteron_m1s_y");
   gProtonSP->SetName("proton_p1s_y");
   gDeuteronSP->SetName("deuteron_p1s_y");

//   TF1* fPionM      = new TF1("fPionM","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fPionM->SetParameters(1.,1.,0.5);
//   TF1* fProtonM    = new TF1("fProtonM","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fProtonM->SetParameters(1.,1.,0.5);
//   TF1* fDeuteronM  = new TF1("fDeuteronM","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fDeuteronM->SetParameters(1.,1.,0.5);
//   TF1* fPionSM     = new TF1("fPionSM","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fPionSM->SetParameters(1.,1.,0.5);
//   TF1* fProtonSM   = new TF1("fProtonSM","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fProtonSM->SetParameters(1.,1.,0.5);
//   TF1* fDeuteronSM = new TF1("fDeuteronSM","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fDeuteronSM->SetParameters(1.,1.,0.5);
//   TF1* fPionSP     = new TF1("fPionSP","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fPionSP->SetParameters(1.,1.,0.5);
//   TF1* fProtonSP   = new TF1("fProtonSP","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fProtonSP->SetParameters(1.,1.,0.5);
//   TF1* fDeuteronSP = new TF1("fDeuteronSP","[0]+[1]/TMath::Power(x,[2])",100.,600.);
//   fDeuteronSP->SetParameters(1.,1.,0.5);

   TF1* fProtonM    = new TF1("fProtonM",   "pol5",0.,500.);
   TF1* fDeuteronM  = new TF1("fDeuteronM", "pol5",0.,500.);
   TF1* fProtonSM   = new TF1("fProtonSM",  "pol5",0.,500.);
   TF1* fDeuteronSM = new TF1("fDeuteronSM","pol5",0.,500.);
   TF1* fProtonSP   = new TF1("fProtonSP",  "pol5",0.,500.);
   TF1* fDeuteronSP = new TF1("fDeuteronSP","pol5",0.,500.);

   fProtonM->SetParameters(1.,1.,1.,1.,1.,1.);
   fProtonSM->SetParameters(1.,1.,1.,1.,1.,1.);
   fProtonSP->SetParameters(1.,1.,1.,1.,1.,1.);

   fProtonM->SetLineColor(2);
   fDeuteronM->SetLineColor(4);
   fProtonSM->SetLineColor(2);
   fDeuteronSM->SetLineColor(4);
   fProtonSP->SetLineColor(2);
   fDeuteronSP->SetLineColor(4);
   fProtonSM->SetLineStyle(2);
   fDeuteronSM->SetLineStyle(2);
   fProtonSP->SetLineStyle(2);
   fDeuteronSP->SetLineStyle(2);

   Double_t fLeft = 0.;
   Double_t fRight = 500.;
   Double_t par[6];

   gProtonM->Fit(fProtonM,"+M0Q","",fLeft,fRight);
   gProtonM->Write();
   fProtonM->SetName("fit_proton_mean_y");
   fProtonM->Write();

   fProtonM->GetParameters(&par[0]);
   fProtonSM->SetParameters(par);
   fProtonSP->SetParameters(par);

   gProtonSM->Fit(fProtonSM,"+M0Q","",fLeft,fRight);
   gProtonSM->Write();
   fProtonSM->SetName("fit_proton_m1s_y");
   fProtonSM->Write();

   gProtonSP->Fit(fProtonSP,"+M0Q","",fLeft,fRight);
   gProtonSP->Write();
   fProtonSP->SetName("fit_proton_1s_y");
   fProtonSP->Write();

   gDeuteronM->Fit(fDeuteronM,"+M0Q","",fLeft,fRight);
   gDeuteronM->Write();
   fDeuteronM->SetName("fit_deuteron_mean_y");
   fDeuteronM->Write();

   fDeuteronM->GetParameters(&par[0]);
   fDeuteronSM->SetParameters(par);
   fDeuteronSP->SetParameters(par);

   gDeuteronSM->Fit(fDeuteronSM,"+M0Q","",fLeft,fRight);
   gDeuteronSM->Write();
   fDeuteronSM->SetName("fit_deuteron_m1s_y");
   fDeuteronSM->Write();

   gDeuteronSP->Fit(fDeuteronSP,"+M0Q","",fLeft,fRight);
   gDeuteronSP->Write();
   fDeuteronSP->SetName("fit_deuteron_1s_y");
   fDeuteronSP->Write();

   TCanvas* d = new TCanvas("dE_TOF_y","dE_TOF_y",800,600);
   d->cd();

   gProtonM->GetXaxis()->SetRangeUser(-50.,50.);
   gProtonM->GetYaxis()->SetRangeUser(0.,500.);
   gProtonM->Draw("AP");
//   gDeuteronM->Draw("Psame");
   gProtonSP->Draw("Psame");
//   gDeuteronSP->Draw("Psame");
   gProtonSM->Draw("Psame");
//   gDeuteronSM->Draw("Psame");

   fProtonM->Draw("same");
//   fDeuteronM->Draw("same");
   fProtonSP->Draw("same");
//   fDeuteronSP->Draw("same");
   fProtonSM->Draw("same");
//   fDeuteronSM->Draw("same");

   d->Write();
   delete d;

   TCanvas* e = new TCanvas("dE_TOF_2D_y","dE_TOF_2D_y",800,600);
   e->cd();

   h->Draw("col");
   fProtonM->SetLineColor(0);
   fProtonSM->SetLineColor(0);
   fProtonSP->SetLineColor(0);
   fProtonM->Draw("same");
   fProtonSM->Draw("same");
   fProtonSP->Draw("same");
   e->Write();

   delete e;

   delete fProtonM;
   delete fDeuteronM;
   delete fProtonSP;
   delete fDeuteronSP;
   delete fProtonSM;
   delete fDeuteronSM;

   delete gProtonM;
   delete gDeuteronM;
   delete gProtonSP;
   delete gDeuteronSP;
   delete gProtonSM;
   delete gDeuteronSM;

   fin.close();

   return;
}

//______________________________________________________________________________

//______________________________________________________________________________
