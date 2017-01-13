/*************************************************************************
 * Author: Manuel Dieterle, 2012
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMFit                                                                //
//                                                                      //
// Class for collecting any kind of TObjects                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TMFit.h"

ClassImp(TMFit)

//______________________________________________________________________________
TMFit::TMFit(const Char_t* szName, const Char_t* szTitle)
    : TNamed(szName, szTitle)
{
    // init fit collection
    fc = new TMFitCollection(szName,szTitle);

    // TMFitCollection config vals
    fNBG = 0;
    fNCTbins = 0;
    fNEbins = 0;
    fRep = 0;

    // TMFit vals
    cData = 0;
    cMCsig = 0;
    cMCbg = 0;
    cEff = 0;
    cBGconfig = 0;
    cBGconfig1 = 0;
    cBGconfig2 = 0;
    cBGindices = 0;
    cCutHisto = 0;
    cExFunc = 0;
    cAxis = 0;
    fRebin = 1;
    fFrac = 1.;
    fTargDens = 0.;
    fBrRatio = 0.;
    fIM = 0.;

    // other vals
    kTwoSig = kFALSE;
    kFitE = kFALSE;
    kExFunc = kFALSE; 
    kIntAll = kFALSE;
    kEff = kTRUE;
    fNBGall = 0;
    hNameBGall = 0;
    fIBG = 0;
    fIBGall = 0;

    cParFitMeth = 0;
    cParFitPolCT = 0;
    cParFitPolE = 0;
    cFitOpt = 0;
    fGloMin = -1.;
    fGloMax = -1.;
    fSigMin = -1.;
    fSigMax = -1.;
    fSigHeightMin = 0.;
    fSigHeightMax = 1.;
    fGauMin = -1.;
    fGauMax = -1.;
    fGauHeightMin = 0.;
    fGauHeightMax = 1.;
    fGauMeanMin = -1.;
    fGauMeanMax = -1.;
    fGauMeanTol = 1.;
    fGauSigmaMin = 0.;
    fGauSigmaMax = 100.;
    fSigma = 3.;
    fMin = 0.;
    fMax = 3000.;
    fFitMin = 0;
    fFitMax = 0;
    cParEPolMean   = 0;
    cParEPolSigma  = 0;
    cParEPolLo     = 0;  
    cParEPolUp     = 0;
    cParCTPolMean  = 0;
    cParCTPolSigma = 0;
    cParCTPolLo    = 0;
    cParCTPolUp    = 0;

    nIt = 0;
    nPadX = 1;
    nPadY = 1;
    ndc = 510;
    xMin = -9999.;
    xMax = 9999.;
}
//______________________________________________________________________________
TMFit::~TMFit()
{
    // Destructor.
    if (fc) delete fc;
}

//______________________________________________________________________________
void TMFit::DrawCuts()
{
   fc->GetCutsFromGaussian();
   DrawCutRegions();
}

//______________________________________________________________________________
void TMFit::DrawFit(Int_t nPadX, Int_t nPadY, Int_t ndc)
{
   Char_t szName[256];

   Int_t iNPad = nPadX*nPadY;
   Int_t nCTbin = fc->GetNumCosThetaBins();
   Int_t nTCbin = fc->GetNumEnergyBins();

   for (Int_t i = 0; i < nCTbin; i++)
   {
      TCanvas* c = 0;
      TPad* p[nTCbin];

      for (Int_t j = 0; j < nTCbin; j++)
      {
         if (j%iNPad == 0)
         {
            if (j+iNPad>nTCbin)
            {
               sprintf(szName, "%s_CT_%i_E_%i_to_%i",     fc->GetName(), i+1, j+1, nTCbin);
               Info("TMFit::DrawFit","Created Pads for CTbin %i TCbin %i to %i",i+1,j+1,nTCbin);
            }
            else
            {
               sprintf(szName, "%s_CT_%i_E_%i_to_%i",     fc->GetName(), i+1, j+1, j+iNPad);
               Info("TMFit::DrawFit","Created Pads for CTbin %i TCbin %i to %i",i+1,j+1,j+iNPad);
            }

            c = new TCanvas(szName,szName,1000,1000);
         }

         c->cd();

         sprintf(szName,"%s_CT_%i_TC_%i",fc->GetName(),i+1,j+1);

         p[j] = new TPad(szName, szName, 0.15+(j%nPadX)*0.8/nPadX, 0.9125-((j%iNPad)/nPadX+1)*0.8/nPadY, 0.15+(j%nPadX+1)*0.8/nPadX, 0.9125-((j%iNPad)/nPadX)*0.8/nPadY);

         p[j]->SetLeftMargin(0);
         p[j]->SetRightMargin(0);
         p[j]->SetTopMargin(0);
         p[j]->SetBottomMargin(0);

         p[j]->Draw();
         p[j]->cd();

         fc->Draw(i,j,xMin,xMax,0,0);

         if (j%iNPad == iNPad-1 || j == nTCbin-1)
         {
            c->cd();

            TGaxis* aX = 0;

            for (Int_t l = 0; l < nPadX; l++)
            {
               aX = new TGaxis(0.15+l*0.8/nPadX,
                                       0.9125-((j%iNPad)/nPadX+1)*0.8/nPadY,
                                       0.15+(l+1)*0.8/nPadX,
                                       0.9125-((j%iNPad)/nPadX+1)*0.8/nPadY,
                                       xMin,
                                       xMax,
                                       ndc,"-");
               aX->SetLabelSize(.03);
               aX->SetLabelOffset(-.03);
               aX->SetLabelFont(42);
               aX->Draw();
            }

            TLine aLine;
            aLine.DrawLineNDC(0.95,0.9125-((j%iNPad)/nPadX+1)*0.8/nPadY,0.95,0.9125);

            TLatex title;
            title.SetTextAlign(22);
            title.SetTextSize(0.05);
            title.DrawLatex(0.15+0.4,0.5*(0.1125-aX->GetLabelSize())+(nPadY-((j%iNPad)/nPadX+1))*0.8/nPadY,Form("#font[41]{%s}",cAxis));
            sprintf(szName,"#font[41]{%2.1f < cos(#theta) < %2.1f}",fc->GetLowerCosTheta(i),fc->GetUpperCosTheta(i));
            title.DrawLatex(0.15+0.4,0.95,szName);
            title.SetTextSize(0.05);
            title.SetTextAngle(90);
            title.DrawLatex(0.5*0.15,0.9125-0.5*((j%iNPad)/nPadX+1)/nPadY*0.8,"#font[41]{Counts [a.u.]}");

            c->Write();
         }
      }
   }

   Info("TMFit::DrawFit","Fits drawn");

   return;
}

//______________________________________________________________________________
void TMFit::DrawCutRegions()
{
   for (Int_t i = 0; i < fc->GetNumCosThetaBins(); i++)
   {
      TCanvas* c = new TCanvas(Form("Data_CT_%i",i),Form("Data_CT_%i",i),1000,1000);
      c->cd();
      fc->DrawInputE(i,"data",cAxis);
      fc->DrawCuts(i,"E");

      c->Update();
      c->Write();
      delete c;

      TCanvas* d = new TCanvas(Form("Signal_CT_%i",i),Form("Signal_CT_%i",i),1000,1000);
      d->cd();
      fc->DrawInputE(i,"signal",cAxis);
      fc->DrawCuts(i,"E");

      d->Update();
      d->Write();
      delete d;
   }

   Info("TMFit::DrawCutRegions","Cut regions (E) drawn");

   for (Int_t i = 0; i < fc->GetNumEnergyBins(); i++)
   {
      TCanvas* c = new TCanvas(Form("Data_E_%i",i),Form("Data_E_%i",i),1000,1000);
      c->cd();
      fc->DrawInputCT(i,"data",cAxis);
      fc->DrawCuts(i,"CT");

      c->Update();
      c->Write();
      delete c;

      TCanvas* d = new TCanvas(Form("Signal_E_%i",i),Form("Signal_E_%i",i),1000,1000);
      d->cd();
      fc->DrawInputCT(i,"signal",cAxis);
      fc->DrawCuts(i,"CT");

      d->Update();
      d->Write();
      delete d;
   }

   Info("TMFit::DrawCutRegions","Cut regions (CT) drawn");


   return;
}


//______________________________________________________________________________
void TMFit::DrawParameters()
{
   for (Int_t i = 0; i < fc->GetNumCosThetaBins(); i++)
   {
      for (Int_t j = 0; j < fc->GetNumBackground()+1; j++)
      {
         TCanvas* c = new TCanvas(Form("ParE_%i_%i",i,j),Form("ParE_%i_%i",i,j),10,10,1200,1200);

         fc->DrawParametersE(i,j);

         c->Write();
         delete c;
      }
   }
   Info("TMFit::DrawParameters","Drawn all energy dependent fit parameters");

   for (Int_t i = 0; i < fc->GetNumEnergyBins(); i++)
   {
      for (Int_t j = 0; j < fc->GetNumBackground()+1; j++)
      {
         TCanvas* c = new TCanvas(Form("ParCT_%i_%i",i,j),Form("ParCT_%i_%i",i,j),10,10,1200,1200);

         fc->DrawParametersCT(i,j);

         c->Write();
         delete c;
      }
   }
   Info("TMFit::DrawParameters","Drawn all cos(theta) dependent fit parameters");

   return;
}

//______________________________________________________________________________
void TMFit::DrawChiSquare()
{
   for (Int_t j = 0; j < fc->GetNumCosThetaBins(); j++)
   {
      TGraphErrors* g = new TGraphErrors(fNEbins);
      for (Int_t k = 0; k < fc->GetNumEnergyBins(); k++)
      {
         Double_t c = fc->GetChiSquare(j,k);
         if (c==TMath::Infinity() || TMath::IsNaN(c)) c = 0.;
         g->SetPoint(k,fc->GetEnergy(k),c);
         g->SetPointError(k,0.,0.);
      }
      fc->SetChiSquareGraphE(j,g);
   }
   for (Int_t j = 0; j < fc->GetNumEnergyBins(); j++)
   {
      TGraphErrors* g = new TGraphErrors(fNCTbins);
      for (Int_t k = 0; k < fc->GetNumCosThetaBins(); k++)
      {
         Double_t c = fc->GetChiSquare(k,j);
         if (c==TMath::Infinity() || TMath::IsNaN(c)) c = 0.;
         g->SetPoint(k,fc->GetLowerCosTheta(k)+0.5*(fc->GetUpperCosTheta(k)-fc->GetLowerCosTheta(k)),c);
         g->SetPointError(k,0.,0.);
      }
      fc->SetChiSquareGraphCT(j,g);
   }

   for (Int_t i = 0; i < fc->GetNumCosThetaBins(); i++)
   {
      TCanvas* c = new TCanvas(Form("ChiE_%i",i),Form("ChiE_%i",i),10,10,1200,1200);

      fc->DrawChiSquareE(i);

      c->Write();
      delete c;
   }
   Info("TMFit::DrawChiSquare","Drawn all energy dependent chi square distributions");

   for (Int_t i = 0; i < fc->GetNumEnergyBins(); i++)
   {
      TCanvas* c = new TCanvas(Form("ChiCT_%i",i),Form("ChiCT_%i",i),10,10,1200,1200);

      fc->DrawChiSquareCT(i);

      c->Write();
      delete c;
   }
   Info("TMFit::DrawChiSquare","Drawn all cos(theta) dependent chi square distributions");

   return;
}

//______________________________________________________________________________
void TMFit::SaveAll()
{
    // init output file
    TFile* fOut = new TFile(Form("%s.root",fc->GetName()),"recreate");
    fOut->cd();
    DrawFit(3,3,504);
    fOut->cd();
    DrawCuts();
    fOut->cd();
    DrawChiSquare();
    fOut->cd();
    DrawParameters();

    fc->CalculateExcitationFunctions(kEff);
    fOut->cd();
    for (Int_t i = 0; i < fNCTbins; i++)
    {
       TH2F* h = (TH2F*)fc->GetExcitationFunction(i);
       h->Write();
    }

    TOEnergyThetaData* sb = (TOEnergyThetaData*)fc->GetSBObject();
    sb->Write();

    TOEnergyThetaCut* fe = (TOEnergyThetaCut*)fc->GetCutObject();
    fe->Write();

    if (fc)
    {
       fc->Write();
       Info("TMFit::SaveAll","Saved Fit Collection %s to file %s.root",fc->GetName(),fc->GetTitle());
    }

    fOut->Close();

    return;

}

//______________________________________________________________________________
void TMFit::CleanGraph(TGraphErrors* g)
{
   Int_t k = g->GetN();
   Int_t n = 0;
   Bool_t kDirty = kTRUE;

   while (kDirty)
   {
      Int_t j = -1;
      for (Int_t i = 0; i < g->GetN(); i++)
      {
         Double_t x = 0.;
         Double_t y = 0.;
         g->GetPoint(i,x,y);

         if (y==0.0)
         {
            j = i;
            n++;
            break;
         }
      }
      if (j >= 0)
         g->RemovePoint(j);
      else
         kDirty = kFALSE;
   }

   Info("TMFit::CleanGraph","Removed %i/%i points of %s %s",n,k,g->Class_Name(),g->GetName());

   return;
}

//______________________________________________________________________________
void TMFit::FitParameters()
{
   for (Int_t i = 0; i < fNCTbins; i++)
   {
      for (Int_t j = 0; j < fNBG+1; j++)
      {
         TGraphErrors* g = new TGraphErrors(fNEbins);

         for (Int_t k = 0; k < fNEbins; k++)
         {
            g->SetPoint(k,fc->GetEnergy(k),fc->GetFitParameter(i,k,j));
            g->SetPointError(k,0.,fc->GetFitParError(i,k,j));
         }
         CleanGraph(g);

         TF1* f = new TF1(Form("fE_%i_%i",i,j),cParFitPolE,fc->GetEnergy(0),fc->GetEnergy(fNEbins-1));
         g->Fit(f,"+RW0MQ");//,"",1000,1500);

         fc->SetEParameterGraph(i,j,g);
         fc->SetEParameterFit(i,j,f);
      }
   }

   for (Int_t i = 0; i < fNEbins; i++)
   {
      for (Int_t j = 0; j < fNBG+1; j++)
      {
         TGraphErrors* g = new TGraphErrors(fNCTbins);

         for (Int_t k = 0; k < fNCTbins; k++)
         {
            g->SetPoint(k,fc->GetLowerCosTheta(k)+0.5*(fc->GetUpperCosTheta(k)-fc->GetLowerCosTheta(k)),
                          fc->GetFitParameter(k,i,j));
            g->SetPointError(k,0.,fc->GetFitParError(k,i,j));
         }
         CleanGraph(g);

         TF1* f = new TF1(Form("fCT_%i_%i",i,j),cParFitPolCT,-1.,1.);
         g->Fit(f,"+RBQ0M");

         fc->SetCTParameterGraph(i,j,g);
         fc->SetCTParameterFit(i,j,f);
      }
   }

}

//______________________________________________________________________________
void TMFit::Fit()
{
   for (Int_t i = 0; i < fRep; i++)
   {
      if (i==0)
         Iterate(i,kTRUE);//,kFitChi);
      else
         Iterate(i,kFALSE);//,kFitChi);
      FitParameters();
   }

   return;
}

//______________________________________________________________________________
void TMFit::Iterate(Int_t n, Bool_t kInit)
{
   for (Int_t i = 0; i < fNCTbins; i++)
   {
      for (Int_t j = 0; j < fNEbins; j++)
      {
         TH1D* hDat = (TH1D*)fc->GetDataObject(i,j);

         TH1** hFit;
         hFit = new TH1*[fNBG+1];
         hFit[0] = (TH1D*)fc->GetSignalObject(i,j);
         
         for (Int_t k = 0; k < fNBG; k++)
            hFit[k+1] = (TH1D*)fc->GetBackgroundObject(i,j,k);
         
         TF1* fFit = 0;
         
         if (kInit)
         {
            // fit only with signal first
            TMFunctor::InitFunctor(1,hFit);
            TF1* fFitTmp = new TF1(Form("fitTmp_%i_%i",i,j),TMFunctor::MyFitFuncMC,fGloMin,fGloMax,1);
            Double_t dMax = 1.;
            if (hFit[0]->GetMaximum()>0.) dMax = hFit[0]->GetMaximum();
            fFitTmp->SetParameter(0,hDat->GetMaximum()/dMax);
            hDat->Fit(fFitTmp,cFitOpt,"",fSigMin,fSigMax);
         
            // fit also with BG contributions now
            TMFunctor::InitFunctor(fNBG+1,hFit);
            fFit = new TF1(Form("fit_%i_%i",i,j),TMFunctor::MyFitFuncMC,fGloMin,fGloMax,fNBG+1);
         
            fFit->SetParameter(0,fFitTmp->GetParameter(0));
            fFit->SetParLimits(0,fSigHeightMin*fFitTmp->GetParameter(0),fSigHeightMax*fFitTmp->GetParameter(0));
            for (Int_t k = 0; k < fNBG; k++)
            {
               Double_t bMax = 1.;
               if (hFit[k+1]->GetMaximum()>0.) bMax = hFit[k+1]->GetMaximum();
               fFit->SetParameter(k+1,hDat->GetMaximum()/bMax);
               fFit->SetParLimits(k+1,1.e-6,1.5*hDat->GetMaximum()/bMax);
            }
         
            hDat->Fit(fFit,cFitOpt);
         }
         else
         {
            // fit with parameter initialization from fit
            TMFunctor::InitFunctor(fNBG+1,hFit);
            fFit = new TF1(Form("fit_%i_%i",i,j),TMFunctor::MyFitFuncMC,fGloMin,fGloMax,fNBG+1);
            for (Int_t k = 0; k < fNBG+1; k++)
            {
               Double_t par = 0.;
               if (kFitE) par = fc->GetEParameterFit(i,k)->Eval(fc->GetEnergy(j));
               else par = fc->GetCTParameterFit(j,k)->Eval(fc->GetLowerCosTheta(i)+0.5*(fc->GetUpperCosTheta(i)-fc->GetLowerCosTheta(i)));
               if (par < 0.0) par = TMath::Abs(par);
               fFit->SetParameter(k,par);
               fFit->SetParLimits(k,fParMin*par,fParMax*par);
            }
            hDat->Fit(fFit,cFitOpt);
         }           

         for (Int_t k = 0; k < fNBG+1; k++)
         {
            fc->SetFitParameter(i,j,k,fFit->GetParameter(k));
            fc->SetFitParError(i,j,k,fFit->GetParError(k));
         }
         fc->SetChiSquare(i,j,fFit->GetChisquare()/fFit->GetNDF());
      }
   }

   return;
}
//______________________________________________________________________________
void TMFit::Configure()
{
   // check if using two signals
   if (strcmp(cBGconfig1,"null") && strcmp(cBGconfig2,"null"))
   {
      kTwoSig = kTRUE;
      Info("TMFit::Configure","Using two signals for fit");
   }
   else
      Info("TMFit::Configure","Using one signal for fit");

   // check if using energy dependent parameter fit
   if (!strcmp(cParFitMeth,"E"))
   {
      kFitE = kTRUE;
      Info("TMFit::Configure","Using energy dependent parameter fit");
   }
   else
      Info("TMFit::Configure","Using cos(theta) dependent parameter fit");

   // set parameter iteration to be applied
   if (!nIt) nIt = fRep;
   Info("TMFit::Configure","Will apply iteration %i/%i to the spectra",nIt,fRep);

   // check if ex funcs should be calculated
   if (strcmp(cEff,"null") && fIM != 0. && fTargDens != 0. && fBrRatio != 0.0)
   {
      kExFunc = kTRUE;
      Info("TMFit::Configure","Excitation function calculation enabled");
   }
   else
   {
      if (strcmp(cEff,"null")) Error("TMFit::Configure","Efficiency file not specified");
      if (fIM == 0.) Error("TMFit::Configure","Invariant mass value not specified");
      if (fTargDens == 0.) Error("TMFit::Configure","Target density not specified");
      if (fBrRatio == 0.) Error("TMFit::Configure","Branching ratio not specified");
      Info("TMFit::Configure","Excitation function calculation disabled");
   }

   // init options for gaussian fit
   fc->SetGaussianOpt(fSigma, fGauMin,fGauMax,fGauMeanMin,fGauMeanMax, fGauMeanTol,
                                      fGauSigmaMin,fGauSigmaMax,
                                      fGauHeightMin,fGauHeightMax);
   Info("TMFit::Configure","Initialized Gaussian Fit Options");

   // set background contributions
   SetBackgroundContributions();

   // set collection
   SetCollection();
}

//______________________________________________________________________________
void TMFit::SetCollection()
{
    // get input files
    TFile* fData = new TFile(cData);
    TFile* fMC = 0;
    TFile* fMC1 = 0;
    TFile* fMC2 = 0;

    if (kTwoSig)
    {
       fMC1 = new TFile(cMCsig1);
       fMC2 = new TFile(cMCsig2);
    }
    else
       fMC = new TFile(cMCsig);   

    TFile** fMCbg;
    fMCbg = new TFile*[fNBG];
    for (Int_t i = 0; i < fNBG; i++)
       fMCbg[i] = new TFile(Form("%s_%i.root",cMCbg,fIBG[i]));

    TFile* fEff = new TFile(cEff);
    TOEnergyThetaData* eff = (TOEnergyThetaData*)fEff->Get("eff_smooth")->Clone("eff");
    fc->SetExFuncCalcOpt(eff,fBrRatio,fTargDens,fIM,cExFunc,kIntAll);
    Info("TMFit::SetCollection","Set ExFunc options, using '%s' as efficiency",fEff->GetName());

    // set dimensions for collection
    TH2F* hTmp = (TH2F*)fData->Get(Form("%s_0",cCutHisto))->Clone("hTmp");
    fNEbins = hTmp->GetNbinsX();
    fc->SetDimensions(fNCTbins, fNEbins, fNBG, fRep, nIt);
    fc->InitObjects();
    fc->SetFitLowerBoundary(fMin,fFitMin);
    fc->SetFitUpperBoundary(fMax,fFitMax);
    fc->SetFitPolMeanE(cParEPolMean);
    fc->SetFitPolSigmaE(cParEPolSigma);
    fc->SetFitPolLoE(cParEPolLo);
    fc->SetFitPolUpE(cParEPolUp);
    fc->SetFitPolMeanCT(cParCTPolMean);
    fc->SetFitPolSigmaCT(cParCTPolSigma);
    fc->SetFitPolLoCT(cParCTPolLo);
    fc->SetFitPolUpCT(cParCTPolUp);

    // set energies
    for (Int_t i = 0; i < hTmp->GetNbinsX(); i++)
       fc->SetEnergy(i,hTmp->GetXaxis()->GetBinCenter(i+1));

    // set cos theta
    for (Int_t i = 0; i < fNCTbins; i++)
    {
       fc->SetLowerCosTheta(i,-1.+i*2./fNCTbins);
       fc->SetUpperCosTheta(i,-1.+(i+1)*2./fNCTbins);
    }

    // set input and fit histograms
    for (Int_t i = 0; i < fNCTbins; i++)
    {
       TH2F* hDatIn = (TH2F*)fData->Get(Form("%s_%i",cCutHisto,i))->Clone(Form("hDataIn_%i",i));
//       hDatIn->Rebin2D(1,fRebin);
       fc->AddInputDataObject(i,hDatIn);

       TH2F* hSigIn = 0;
       if (kTwoSig)
       {
	  hSigIn = (TH2F*)fMC1->Get(Form("%s_%i",cCutHisto,i))->Clone(Form("hSigIn_%i",i));
	  hSigIn->Add((TH2F*)fMC2->Get(Form("%s_%i",cCutHisto,i)), fFrac);
       }
       else
	  hSigIn = (TH2F*)fMC->Get(Form("%s_%i",cCutHisto,i))->Clone(Form("hSigIn_%i",i));
//       hSigIn->Rebin2D(1,fRebin);
       fc->AddInputSignalObject(i,hSigIn);

       TH2F* hBGIn[fNBG];
       for (Int_t j = 0; j < fNBG; j++)
       {
          hBGIn[j] = (TH2F*)fMCbg[j]->Get(Form("%s_%i",cCutHisto,i))->Clone(Form("hBGIn_%i_%i",i,j));
//	  hBGIn[j]->Rebin2D(1,fRebin);
          fc->AddInputBackgroundObject(i,j,hBGIn[j]);
       }
       for (Int_t k = 0; k < hTmp->GetNbinsX(); k++)
       {
	  TH1D* hDat = (TH1D*)hDatIn->ProjectionY(Form("hDat_%i_%i",i,k),k+1,k+1);
          hDat->Rebin(fRebin);
	  if (!strcmp(cAxis,"null"))
	     hDat->GetXaxis()->SetTitle("axis title not set");
	  else
	     hDat->GetXaxis()->SetTitle(cAxis);
	  hDat->GetXaxis()->CenterTitle();
	  hDat->GetYaxis()->SetTitle("#font[42]{Counts [a.u.]}");
	  hDat->GetYaxis()->CenterTitle();
          fc->AddDataObject(i,k,hDat);
	  TH1D* hSig = (TH1D*) hSigIn ->ProjectionY(Form("hSig_%i_%i",i,k),k+1,k+1);
          hSig->Rebin(fRebin);
          fc->AddSignalObject(i,k,hSig);
	  for (Int_t j = 0; j < fNBG; j++)
          {
	     TH1D* hBG = (TH1D*)hBGIn[j]->ProjectionY(Form("hBG_%i_%i_%i",i,k,j),k+1,k+1);
             hBG->Rebin(fRebin);
             fc->AddBackgroundObject(i,k,j,hBG);
          }
       }
    }

    Info("TMFit::SetCollection","Set Collection!");

    return;
}


//______________________________________________________________________________
void TMFit::SetBackgroundContributions()
{
   // count number of BG contributions
//   fNBGall = 0;

   // if using two signals read both files
   if (kTwoSig)
   {
      ifstream f(cBGconfig1);
      string line1;
      for (Int_t i = 0; getline(f, line1); i++)
	 fNBGall++;
      f.clear();
      f.seekg(0);
      f.close();

      ifstream g(cBGconfig2);
      string line2;
      for (Int_t i = 0; getline(g, line2); i++)
	 fNBGall++;
      g.clear();
      g.seekg(0);
      g.close();
   }
   // otherwise only one
   else
   {
      ifstream f(cBGconfig);
      string line;
      for (Int_t i = 0; getline(f, line); i++)
	 fNBGall++;
      f.clear();
      f.seekg(0);
      f.close();
   }

   // count bg contributions to use
   fNBG = 0;
   if (cBGindices)
   {
      for (UInt_t i = 0; i < strlen(cBGindices); i++)
	 if (cBGindices[i] == ',') fNBG++;
      fNBG++;
   }

   // get bg contributions to use
   fIBG = new Int_t[fNBG];
   Int_t pos = 0;
   if (cBGindices)
   {
      Char_t* tmp;
      tmp = strtok(cBGindices, ",");
      while (tmp != NULL)
      {
	 fIBG[pos] = atoi(tmp);
	 pos++;
	 tmp = strtok(NULL, ",");
      }
   }

   // get background descriptions
   fIBGall = new Int_t[fNBGall];
   hNameBGall = new Char_t*[fNBGall];
  
   if (kTwoSig)
   {
      Int_t n = 0;
      ifstream f(cBGconfig1);
      string line1;
      for (Int_t i = 0; getline(f, line1); i++)
      {
	 Char_t tmp[256];
	 sprintf(tmp,"%s",line1.c_str());
	 fIBGall[n] = atoi(strtok(tmp, " "));
	 Char_t* ttmp = strtok(NULL, " ");
	 hNameBGall[n] = new Char_t[128];
	 sprintf(hNameBGall[n],"%s",ttmp);
	 n++;
      }
      f.clear();
      f.seekg(0);
      f.close();
      ifstream g(cBGconfig2);
      string line2;
      for (Int_t j = 0; getline(g, line2); j++)
      {
	 Char_t tmp[256];
	 sprintf(tmp,"%s",line2.c_str());
	 fIBGall[n] = atoi(strtok(tmp, " "));
	 Char_t* ttmp = strtok(NULL, " ");
	 hNameBGall[n] = new Char_t[128];
	 sprintf(hNameBGall[n],"%s",ttmp);
	 n++;
      }
      g.clear();
      g.seekg(0);
      g.close();
   }
   else
   {
      ifstream f(cBGconfig);
      string line;
      for (Int_t i = 0; getline(f, line); i++)
      {
	 Char_t tmp[256];
	 sprintf(tmp,"%s",line.c_str());
	 fIBGall[i] = atoi(strtok(tmp, " "));
	 Char_t* ttmp = strtok(NULL, " ");
	 hNameBGall[i] = new Char_t[128];
	 sprintf(hNameBGall[i],"%s",ttmp);
      }
      f.clear();
      f.seekg(0);
      f.close();
   }

   return;
}


//______________________________________________________________________________
void TMFit::Print()
{
   Printf("========================================================");
   Printf("|| TMFit Settings:");
   Printf("========================================================");
   Printf("|| Number of CT bins:                   %i", fNCTbins);
   Printf("|| Number of E bins:                    %i", fNEbins);
   Printf("|| Number of Background Contributions:  %i(%i)", fNBG, fNBGall);
   Printf("|| Number of Fit Repetitions:           %i", fRep);
   Printf("|| Input Data file:                     %s", cData);
   if (kTwoSig)
   {
      Printf("|| Input MC signal1 file:               %s", cMCsig1);
      Printf("|| Input MC signal2 file:               %s", cMCsig2);
   }
   else
      Printf("|| Input MC signal file:                %s", cMCsig);
   Printf("|| Input MC background files prefix:    %s", cMCbg);
   if (kTwoSig)
   {
      Printf("|| Background1 Config file:             %s", cBGconfig1);
      Printf("|| Background2 Config file: (%2.1f%%)      %s", fFrac, cBGconfig2);
   }
   else
      Printf("|| Background Config file:              %s", cBGconfig);
   Printf("|| Background Contributions:");
   for (Int_t i = 0; i < fNBGall; i++)
   {
      Bool_t kFound = kFALSE;
      for (Int_t j = 0; j < fNBG; j++)
	 if (fIBGall[i]==fIBG[j])
	    kFound = kTRUE;
      if (kFound)
	 Printf("||                                      %i: %s",i , hNameBGall[i]);
      else
	 Printf("||                                      %i: %s (omitted)",i , hNameBGall[i]);
   }

   Printf("|| Parameter Fit Method:                %s", cParFitMeth);
   Printf("|| CT parameter polynomial:             %s", cParFitPolCT);
   Printf("|| E parameter polynomial:              %s", cParFitPolE);
   Printf("|| Fit parameter variation range        [%2.1f,%2.1f]",fParMin, fParMax);
   Printf("|| Fit Option:                          %s", cFitOpt);
   Printf("|| Global Fit Range:                    [%2.1f,%2.1f]",fGloMin, fGloMax);
   Printf("|| Signal Fit Range:                    [%2.1f,%2.1f]",fSigMin, fSigMax);
   Printf("|| Signal Height Range:                 [%2.1f,%2.1f]",fSigHeightMin, fSigHeightMax);
   Printf("|| Gaussian Fit Range:                  [%2.1f,%2.1f]",fGauMin, fGauMax);
   Printf("|| Gaussian Height Range:               [%2.1f,%2.1f]",fGauHeightMin, fGauHeightMax);
   Printf("|| Gaussian Mean Range:                 [%2.1f,%2.1f]",fGauMeanMin, fGauMeanMax);
   Printf("|| Gaussian Mean Tolerance:              %3.0f",fGauMeanTol);
   Printf("|| Gaussian Sigma Range:                [%2.1f,%2.1f]",fGauSigmaMin, fGauSigmaMax);
   Printf("|| Fit Sigma:                           %2.1f",fSigma);
   Printf("|| Cut Range:                           [%2.1f,%2.1f]",fMin, fMax);

   printf("||      Fit Lower Boundaries: [");
   for (Int_t i = 0; i < fNCTbins-1; i++)
      printf("%4.0f, ",fc->GetFitLowerBoundary(i));
   printf("%4.0f]\n",fc->GetFitLowerBoundary(fNCTbins-1));
   printf("||      Fit Upper Boundaries: [");
   for (Int_t i = 0; i < fNCTbins-1; i++)
      printf("%4.0f, ",fc->GetFitUpperBoundary(i));
   printf("%4.0f]\n",fc->GetFitUpperBoundary(fNCTbins-1));

   Printf("|| Cut Fit Range [CT] lo:               %s",fFitMin);
   Printf("|| Cut Fit Range [CT] up:               %s",fFitMax);

   printf("|| E - Fit Pol Mean: [");
   for (Int_t i = 0; i < fNCTbins-1; i++)
      printf("%s, ",fc->GetFitPolMeanE(i));
   printf("%s]\n",fc->GetFitPolMeanE(fNCTbins-1));
   printf("|| E - Fit Pol Sigma: [");
   for (Int_t i = 0; i < fNCTbins-1; i++)
      printf("%s, ",fc->GetFitPolSigmaE(i));
   printf("%s]\n",fc->GetFitPolSigmaE(fNCTbins-1));
   printf("|| E - Fit Pol Lo: [");
   for (Int_t i = 0; i < fNCTbins-1; i++)
      printf("%s, ",fc->GetFitPolLoE(i));
   printf("%s]\n",fc->GetFitPolLoE(fNCTbins-1));
   printf("|| E - Fit Pol Up: [");
   for (Int_t i = 0; i < fNCTbins-1; i++)
      printf("%s, ",fc->GetFitPolUpE(i));
   printf("%s]\n",fc->GetFitPolUpE(fNCTbins-1));

   Printf("|| CT - Parameter Mean Fit Polynomial:  %s", cParCTPolMean);
   Printf("|| CT - Parameter Sigma Fit Polynomial: %s", cParCTPolSigma);
   Printf("|| CT - Parameter Lo Fit Polynomial:    %s", cParCTPolLo);
   Printf("|| CT - Parameter Up Fit Polynomial:    %s", cParCTPolUp);

   Printf("|| Name of Cut histogram:               %s", cCutHisto);
   Printf("|| Rebin Value for cut histogram:       %i", fRebin);
   Printf("|| Axis title:                          %s", cAxis);

   if (kExFunc)
   {
      Printf("|| Name of Ex Func histogram:           %s", cExFunc);
      printf("|| Efficiency file:                     %s", cEff);
      if (!kEff) printf(" (not applied)\n");
      else printf("\n");
      Printf("|| Target density:                      %f", fTargDens);
      Printf("|| Branching ratio:                     %f", fBrRatio);
      Printf("|| Invariant mass value:                %3.0f", fIM);
   }

   Printf("========================================================");


}

//______________________________________________________________________________
