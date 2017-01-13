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


#include "TMFitCollection.h"

ClassImp(TMFitCollection)

//______________________________________________________________________________
TMFitCollection::TMFitCollection(const Char_t* szName, const Char_t* szTitle)
    : TNamed(szName, szTitle)
{
    // Constructor. Allocate space for 'nMaxTaggCh' tagger channels.
 
    // init dimensions
    fNBG = 0;
    fNCTbins = 1;
    fNEbins = 1;
    fRep = 1;

    fGmin = -9999.;
    fGmax = 9999.;
    fGMmin = -9999.;
    fGMmax = 9999.;
    fGMtol = 1.;
    fGSmin = 0.;
    fGSmax = 9999.;
    fGHmin = 0.;
    fGHmax = 999999999.;

    fSigma = 0.;
    fMin = 0.;
    fMax = 9999.;
    fLo = 0;
    fUp = 0;

    szPolMeanCT = 0;
    szPolSigmaCT = 0;
    szPolLoCT = 0;
    szPolUpCT = 0;

    eCut = 0;
    eff = 0;
    fBR = 1.;
    fTD = 1.;
    fIM = 0.;
    sExFunc = 0;
    kAll = kFALSE;

    sb = 0;
}

//______________________________________________________________________________
TMFitCollection::~TMFitCollection()
{
    if (fLo) delete [] fLo;
    if (fUp) delete [] fUp;
}

//______________________________________________________________________________
void TMFitCollection::DrawChiSquareE(Int_t c)
{
   gChiSquareE[c]->SetTitle("");
   gChiSquareE[c]->SetMarkerColor(kBlack);
   gChiSquareE[c]->SetMarkerSize(1.2);
   gChiSquareE[c]->SetMarkerStyle(20);
   gChiSquareE[c]->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
   gChiSquareE[c]->GetYaxis()->SetTitle(Form("Chi2/ndf (%i rep.) [a.u.]",fRep));
   gChiSquareE[c]->GetXaxis()->CenterTitle();
   gChiSquareE[c]->GetYaxis()->CenterTitle();
   gChiSquareE[c]->Draw("ap");

   return;
}

//______________________________________________________________________________
void TMFitCollection::DrawChiSquareCT(Int_t e)
{
   gChiSquareCT[e]->SetTitle("");
   gChiSquareCT[e]->SetMarkerColor(kBlack);
   gChiSquareCT[e]->SetMarkerSize(1.2);
   gChiSquareCT[e]->SetMarkerStyle(20);
   gChiSquareCT[e]->GetXaxis()->SetTitle("cos(#theta) [a.u.]");
   gChiSquareCT[e]->GetYaxis()->SetTitle(Form("Chi2/ndf (%i rep.) [a.u.]",fRep));
   gChiSquareCT[e]->GetXaxis()->CenterTitle();
   gChiSquareCT[e]->GetYaxis()->CenterTitle();
   gChiSquareCT[e]->Draw("ap");

   return;
}

//______________________________________________________________________________
void TMFitCollection::DrawInputE(Int_t c, const Char_t* s, const Char_t* t)
{
   TH2F* h = 0;

   if (!strcmp(s,"data")) h = (TH2F*)oDatIn[c];
   else if (!strcmp(s,"signal")) h = (TH2F*)oSigIn[c];
   else Error("TMFitCollection::DrawInputE","unknown input type '%s'",s);

   h->GetXaxis()->SetTitle("Eg [MeV]");
   h->GetYaxis()->SetTitle(t);
   h->GetXaxis()->CenterTitle();
   h->GetYaxis()->CenterTitle();
   h->Draw("colz");

   return;
}

//______________________________________________________________________________
void TMFitCollection::DrawInputCT(Int_t e, const Char_t* s, const Char_t* t)
{
   TH2F* hTmp = 0;
   if (!strcmp(s,"data")) hTmp = (TH2F*)oDatIn[0];
   else if (!strcmp(s,"signal")) hTmp = (TH2F*)oSigIn[0];
   else Error("TMFitCollection::DrawInputCT","unknown input type '%s'",s);

   TH2F* h = new TH2F(Form("%s_%i",hTmp->GetName(),e),Form("%s_%i",hTmp->GetName(),e),
        fNCTbins,-1,1,
        hTmp->GetNbinsY(),hTmp->GetYaxis()->GetBinLowEdge(1),hTmp->GetYaxis()->GetBinUpEdge(hTmp->GetNbinsY()));
   h->Sumw2();

   TH2F** hct = new TH2F*[fNCTbins];

   for (Int_t i = 0; i < fNCTbins; i++)
   {
      if (!strcmp(s,"data")) hct[i] = (TH2F*)oDatIn[i];
      else if (!strcmp(s,"signal")) hct[i] = (TH2F*)oSigIn[i];
      else Error("TMFitCollection::DrawInputCT","unknown input type '%s'",s);
      for (Int_t j = 0; j < hct[i]->GetNbinsY(); j++)
      {
         h->SetBinContent(i+1,j+1,hct[i]->GetBinContent(e+1,j+1));
         h->SetBinError(i+1,j+1,hct[i]->GetBinError(e+1,j+1));
      }
   }

   h->GetXaxis()->SetTitle("cos(#theta) [rad]");
   h->GetYaxis()->SetTitle(t);
   h->GetXaxis()->CenterTitle();
   h->GetYaxis()->CenterTitle();
   h->Draw("colz");

   return;
}

//______________________________________________________________________________
void TMFitCollection::DrawCuts(Int_t c, const Char_t* s)
{
   if (!strcmp(s,"E"))
   {
      gCutLoE[c]->SetMarkerStyle(20);
      gCutLoE[c]->SetMarkerSize(1.2);
      gCutLoE[c]->SetMarkerColor(kBlack);
      gCutLoE[c]->Draw("lsame");
   
      fCutLoE[c]->SetLineColor(kRed);
      fCutLoE[c]->SetLineWidth(2);
      fCutLoE[c]->Draw("csame");
   
      gCutUpE[c]->SetMarkerStyle(20);
      gCutUpE[c]->SetMarkerSize(1.2);
      gCutUpE[c]->SetMarkerColor(kBlack);
      gCutUpE[c]->Draw("lsame");
   
      fCutUpE[c]->SetLineColor(kRed);
      fCutUpE[c]->SetLineWidth(2);
      fCutUpE[c]->Draw("csame");
   }
   else if (!strcmp(s,"CT"))
   {
      gCutLoCT[c]->SetMarkerStyle(20);
      gCutLoCT[c]->SetMarkerSize(1.2);
      gCutLoCT[c]->SetMarkerColor(kBlack);
      gCutLoCT[c]->Draw("lsame");
   
      fCutLoCT[c]->SetLineColor(kRed);
      fCutLoCT[c]->SetLineWidth(2);
      fCutLoCT[c]->Draw("csame");
   
      gCutUpCT[c]->SetMarkerStyle(20);
      gCutUpCT[c]->SetMarkerSize(1.2);
      gCutUpCT[c]->SetMarkerColor(kBlack);
      gCutUpCT[c]->Draw("lsame");
   
      fCutUpCT[c]->SetLineColor(kRed);
      fCutUpCT[c]->SetLineWidth(2);
      fCutUpCT[c]->Draw("csame");
   }
   else
      Error("TMFit::DrawCuts","unknown type '%s'",s);

   return;
}

//______________________________________________________________________________
void TMFitCollection::GetCutsFromGaussian()
{
   eCut = new TOEnergyThetaCut(fNCTbins,"cut");
   
   for (Int_t i = 0; i < fNCTbins; i++)
   {
      for (Int_t j = 0; j < fNEbins; j++)
      {
         TF1* f = (TF1*)oFit[i][j]->Clone(Form("f_%i_%i",i,j));
         Double_t mean = f->GetParameter(1);
         Double_t mean_err = f->GetParError(1);
         Double_t sigm = f->GetParameter(2);
         Double_t sigm_err = f->GetParError(2);
         gCutMeanE[i]->SetPoint(j,fEnergy[j],mean);
         gCutSigmaE[i]->SetPoint(j,fEnergy[j],sigm);
         gCutLoE[i]->SetPoint(j,fEnergy[j],mean-fSigma*sigm);
         gCutUpE[i]->SetPoint(j,fEnergy[j],mean+fSigma*sigm);
         gCutMeanE[i]->SetPointError(j,0.,mean_err);
         gCutSigmaE[i]->SetPointError(j,0.,sigm_err);
         gCutLoE[i]->SetPointError(j,0.,TMath::Sqrt(mean_err*mean_err+fSigma*fSigma*sigm_err*sigm_err));
         gCutUpE[i]->SetPointError(j,0.,TMath::Sqrt(mean_err*mean_err+fSigma*fSigma*sigm_err*sigm_err));
      }

      gCutMeanE[i]  = CleanGraphE(i,gCutMeanE[i]);
      gCutSigmaE[i] = CleanGraphE(i,gCutSigmaE[i]);
      gCutLoE[i]    = CleanGraphE(i,gCutLoE[i]);
      gCutUpE[i]    = CleanGraphE(i,gCutUpE[i]);

      fCutMeanE[i] = new TF1(Form("fMeanE_%i",i),szPolMeanE[i],fMin,fMax);
      SetPolynomialParameters(fCutMeanE[i],szPolMeanE[i]);
      gCutMeanE[i]->SetName(Form("gMeanE_%i",i));
      gCutMeanE[i]->Fit(fCutMeanE[i],"+WBQ0M","",fLo[i],fUp[i]);
      fCutSigmaE[i] = new TF1(Form("fSigmaE_%i",i),szPolSigmaE[i],fMin,fMax);
      SetPolynomialParameters(fCutSigmaE[i],szPolSigmaE[i]);
      gCutSigmaE[i]->SetName(Form("gSigmaE_%i",i));
      gCutSigmaE[i]->Fit(fCutSigmaE[i],"+WBQ0M","",fLo[i],fUp[i]);

      fCutLoE[i] = new TF1(Form("fLoE_%i",i),szPolLoE[i],fMin,fMax);
      SetPolynomialParameters(fCutLoE[i],szPolLoE[i]);
      gCutLoE[i]->SetName(Form("gLoE_%i",i));
      gCutLoE[i]->Fit(fCutLoE[i],"+WBQ0M","",fLo[i],fUp[i]);
      fCutUpE[i] = new TF1(Form("fUpE_%i",i),szPolUpE[i],fMin,fMax);
      SetPolynomialParameters(fCutUpE[i],szPolUpE[i]);
      gCutUpE[i]->SetName(Form("gUpE_%i",i));
      gCutUpE[i]->Fit(fCutUpE[i],"+WBQ0M","",fLo[i],fUp[i]);

      eCut->SetMeanGraph(i,gCutMeanE[i]);
      eCut->SetSigmaGraph(i,gCutSigmaE[i]);
      eCut->SetLowerGraph(i,gCutLoE[i]);
      eCut->SetUpperGraph(i,gCutUpE[i]);
      eCut->SetMeanFunction(i,fCutMeanE[i]);
      eCut->SetSigmaFunction(i,fCutSigmaE[i]);
      eCut->SetLowerFunction(i,fCutLoE[i]);
      eCut->SetUpperFunction(i,fCutUpE[i]);
   }

   for (Int_t i = 0; i < fNEbins; i++)
   {
      for (Int_t j = 0; j < fNCTbins; j++)
      {
         TF1* f = (TF1*)oFit[j][i]->Clone(Form("f_ct_%i_%i",j,i));
         Double_t mean = f->GetParameter(1);
         Double_t mean_err = f->GetParError(1);
         Double_t sigm = f->GetParameter(2);
         Double_t sigm_err = f->GetParError(2);
         gCutMeanCT[i]->SetPoint(j,0.5*(fCTlo[j]+fCTup[j]),mean);
         gCutSigmaCT[i]->SetPoint(j,0.5*(fCTlo[j]+fCTup[j]),sigm);
         gCutLoCT[i]->SetPoint(j,0.5*(fCTlo[j]+fCTup[j]),mean-fSigma*sigm);
         gCutUpCT[i]->SetPoint(j,0.5*(fCTlo[j]+fCTup[j]),mean+fSigma*sigm);
         gCutMeanCT[i]->SetPointError(j,0.,mean_err);
         gCutSigmaCT[i]->SetPointError(j,0.,sigm_err);
         gCutLoCT[i]->SetPointError(j,0.,TMath::Sqrt(mean_err*mean_err+fSigma*fSigma*sigm_err*sigm_err));
         gCutUpCT[i]->SetPointError(j,0.,TMath::Sqrt(mean_err*mean_err+fSigma*fSigma*sigm_err*sigm_err));
      }

      CleanGraphCT(i,gCutMeanCT[i]);
      CleanGraphCT(i,gCutSigmaCT[i]);
      CleanGraphCT(i,gCutLoCT[i]);
      CleanGraphCT(i,gCutUpCT[i]);

      fCutMeanCT[i] = new TF1(Form("fMeanCT_%i",i),szPolMeanCT,-1,1);
      SetPolynomialParameters(fCutMeanCT[i],szPolMeanCT);
      gCutMeanCT[i]->SetName(Form("gMeanCT_%i",i));
      gCutMeanCT[i]->Fit(fCutMeanCT[i],"+WBQ0M","",-1,1);
      fCutSigmaCT[i] = new TF1(Form("fSigmaCT_%i",i),szPolSigmaCT,-1,1);
      SetPolynomialParameters(fCutSigmaCT[i],szPolSigmaCT);
      gCutSigmaCT[i]->SetName(Form("gSigmaCT_%i",i));
      gCutSigmaCT[i]->Fit(fCutSigmaCT[i],"+WBQ0M","",-1,1);

      fCutLoCT[i] = new TF1(Form("fLoCT_%i",i),szPolLoCT,-1,1);
      SetPolynomialParameters(fCutLoCT[i],szPolLoCT);
      gCutLoCT[i]->SetName(Form("gLoCT_%i",i));
      gCutLoCT[i]->Fit(fCutLoCT[i],"+WBQ0M","",-1,1);
      fCutUpCT[i] = new TF1(Form("fUpCT_%i",i),szPolUpCT,-1,1);
      SetPolynomialParameters(fCutUpCT[i],szPolUpCT);
      gCutUpCT[i]->SetName(Form("gUpCT_%i",i));
      gCutUpCT[i]->Fit(fCutUpCT[i],"+WBQ0M","",-1,1);
   }

   return;
}

//______________________________________________________________________________
void TMFitCollection::SetPolynomialParameters(TF1* f, Char_t* s)
{
   Int_t p = 0;
   sscanf(s,"pol%i",&p);

   for (Int_t i = 0; i < p+1; i++)
      f->SetParameters(i,1.);

   return;
}

//______________________________________________________________________________
void TMFitCollection::SetTotalHistos(Int_t ct, Int_t e)
{
   TH1D* h[fNBG+1];

   h[0] = (TH1D*)oSig[ct][e]->Clone(Form("hTmp_%i_%i",ct,e));
   for (Int_t k = 0; k < fNBG; k++)
      h[k+1] = (TH1D*)oBG[ct][e][k]->Clone(Form("hTmp_%i_%i_%i",ct,e,k+1));

   for (Int_t k = 0; k < fNBG+1; k++)
   {
      if (k==0)
         oTot[ct][e] = (TH1D*)h[0]->Clone(Form("hTot_%i_%i",ct,e));
      else
      {
         oTot[ct][e]->Add(h[k]);
         if (k==1) oBGTot[ct][e] = (TH1D*)h[k]->Clone(Form("hBGTot_%i_%i",ct,e));
         else oBGTot[ct][e]->Add(h[k]);
      }
   }

   Double_t* g = GetMaximumAndPositionInRange(h[0]); 
//   Double_t mTmp = h[0]->GetBinCenter(h[0]->GetMaximumBin());
//   oFit[ct][e]->SetParameters(h[0]->GetMaximum(),fGMmin+0.5*(fGMmax-fGMmin),fGSmin+0.5*(fGSmax-fGSmin));
//Printf("Max:%f | %f < Mean:%f < %f",g[0],g[2],g[1],g[3]);
   oFit[ct][e]->SetParameters(g[0],g[1],fGSmin+0.5*(fGSmax-fGSmin));
   oFit[ct][e]->SetParLimits(0,fGHmin*g[0],fGHmax*g[0]);
//   oFit[ct][e]->SetParLimits(1,fGMmin,fGMmax);
   oFit[ct][e]->SetParLimits(1,g[2],g[3]);
   oFit[ct][e]->SetParLimits(2,fGSmin,fGSmax);
   h[0]->Fit(oFit[ct][e],"+RBQ0M");

   return;   
}
//______________________________________________________________________________
void TMFitCollection::Draw(Int_t ct, Int_t e, Double_t rLo, Double_t rUp, Bool_t kTit, Bool_t kPrint)
{
   oSig[ct][e]->Scale(fPar[ct][e][0]);
   if (kPrint) 
      printf("Parameter: CT%i(%i) E%i(%i) Rep:%i Sig:%f%%",
             ct+1,fNCTbins,e+1,fNEbins,fRep,oSig[ct][e]->Integral()/oDat[ct][e]->Integral()*100.);

   for (Int_t k = 0; k < fNBG; k++)
   {
      oBG[ct][e][k]->Scale(fPar[ct][e][k+1]);
      if (kPrint)
      {
         printf(" BG%i:%f%%",k+1,oBG[ct][e][k]->Integral()/oDat[ct][e]->Integral()*100.);
         if (k==fNBG-1) printf("\n");
      }
   }
 
   SetTotalHistos(ct,e);
 
   TH1* hDat = (TH1*)oDat[ct][e];
   hDat->SetMarkerColor(kBlack);
   hDat->SetMarkerSize(1.1);
   hDat->SetLineColor(kBlack);

   TH1* hSig = (TH1*)oSig[ct][e];
   hSig->SetLineColor(kBlue);
   hSig->SetLineWidth(2);

   TH1* hBG[fNBG];
   for (Int_t i = 0; i < fNBG; i++)
   {
      hBG[i] = (TH1*)oBG[ct][e][i];
      hBG[i]->SetLineColor(kGreen);
      hBG[i]->SetLineWidth(2);
   }

   TH1* hBGTot = (TH1*)oBGTot[ct][e];
   hBGTot->SetLineColor(kMagenta);
   hBGTot->SetLineWidth(2);

   TH1* hTot = (TH1*)oTot[ct][e];
   hTot->SetLineColor(kRed);
   hTot->SetLineWidth(2);

   TF1* fFit = (TF1*)oFit[ct][e];
   fFit->SetLineColor(kBlue);
   fFit->SetLineStyle(2);
   fFit->SetLineWidth(2);

//   TLine* lLo = new TLine(fCutLow[ct]->Eval(fEnergy[e]),0.0,fCutLow[ct]->Eval(fEnergy[e]),1.2*hDat->GetMaximum());
//   lLo->SetLineStyle(2);
//   lLo->SetLineWidth(2);
//   TLine* lUp = new TLine(fCutUpp[ct]->Eval(fEnergy[e]),0.0,fCutUpp[ct]->Eval(fEnergy[e]),1.2*hDat->GetMaximum());
//   lUp->SetLineStyle(2);
//   lUp->SetLineWidth(2);

   if (kTit)
   {
      hDat->SetTitle(Form("Eg = %4.0f MeV | %2.1f < cos(#theta) < %2.1f | Chi2/ndf=%f",fEnergy[e],fCTlo[ct],fCTup[ct],fChiSquare[ct][e]));
   }
   else
   {
      hDat->SetTitle(Form("Eg = %4.0f MeV | Chi2/ndf=%f",fEnergy[e],fChiSquare[ct][e]));
   }

   hDat->SetStats(0);
   hDat->GetYaxis()->SetRangeUser(0.0,1.2*hDat->GetMaximum());
   hDat->Draw("E1");
   if (rLo != 9999. && rUp != 9999.)
      hDat->GetXaxis()->SetRangeUser(rLo,rUp);
   hSig->Draw("Hsame");
   for (Int_t i = 0; i < fNBG; i++)
      hBG[i]->Draw("Hsame");
   hBGTot->Draw("Hsame");
   hTot->Draw("Hsame");
   fFit->Draw("csame");
//   lLo->Draw("csame");
//   lUp->Draw("csame");
}
//______________________________________________________________________________
void TMFitCollection::DrawParametersCT(Int_t e, Int_t bg)
{
   gParCT[e][bg]->SetTitle("");
   gParCT[e][bg]->SetMarkerColor(kBlack);
   gParCT[e][bg]->SetMarkerSize(1.2);
   gParCT[e][bg]->SetMarkerStyle(20);
   gParCT[e][bg]->GetXaxis()->SetTitle("cos(#theta) [a.u.]");
   gParCT[e][bg]->GetYaxis()->SetTitle(Form("Par_%i (%i rep.) [a.u.]",bg,fRep));
   gParCT[e][bg]->GetXaxis()->CenterTitle();
   gParCT[e][bg]->GetYaxis()->CenterTitle();
   gParCT[e][bg]->Draw("ap");

   fFitCT[e][bg]->SetLineColor(kRed);
   fFitCT[e][bg]->SetLineWidth(2);
   fFitCT[e][bg]->Draw("csame");

   return;
}

//______________________________________________________________________________
void TMFitCollection::DrawParametersE(Int_t ct, Int_t bg)
{
   gParE[ct][bg]->SetTitle("");
   gParE[ct][bg]->SetMarkerColor(kBlack);
   gParE[ct][bg]->SetMarkerSize(1.2);
   gParE[ct][bg]->SetMarkerStyle(20);
   gParE[ct][bg]->GetXaxis()->SetTitle("Eg [MeV]");
   gParE[ct][bg]->GetYaxis()->SetTitle(Form("Par_%i (%i rep.) [a.u.]",bg,fRep));
   gParE[ct][bg]->GetXaxis()->CenterTitle();
   gParE[ct][bg]->GetYaxis()->CenterTitle();
   gParE[ct][bg]->Draw("ap");

   fFitE[ct][bg]->SetLineColor(kRed);
   fFitE[ct][bg]->SetLineWidth(2);
   fFitE[ct][bg]->Draw("csame");

   return;
}

//______________________________________________________________________________
Double_t* TMFitCollection::GetMaximumAndPositionInRange(TH1D* h)
{
   Double_t* p = new Double_t[4];
   Double_t v = h->GetBinCenter(h->GetMaximumBin());

   if (v > fGMmin && v < fGMmax)
   {
      p[0] = h->GetBinContent(h->GetMaximumBin());
      p[1] = v;
   }
   else
   {
      Double_t a = 0;
      Double_t m = 0;
      for (Int_t i = h->FindBin(fGMmin); i < h->FindBin(fGMmax)+1; i++)
      {
         Double_t c = h->GetBinContent(i);
         if (c > a)
         {
            a = c;
            m = h->GetBinCenter(i);
         }
      }

      p[0] = a;
      p[1] = m;
   }

   p[2] = p[1]-fGMtol*TMath::Abs(p[1]);
   p[3] = p[1]+fGMtol*TMath::Abs(p[1]);

   return p;
}

//______________________________________________________________________________
void TMFitCollection::SetFitLowerBoundary(Double_t v, Char_t* lo)
{
   fMin = v;

   Char_t* pch;
   pch = strtok (lo,",");

   Int_t n = 0;
   while (pch != NULL)
   {
      fLo[n] = atof(pch);
      n++;
      pch = strtok (NULL, ", ");
   }

   return;
}

//______________________________________________________________________________
void TMFitCollection::SetFitUpperBoundary(Double_t v, Char_t* up)
{
   fMax = v;

   Char_t* pch;
   pch = strtok (up,",");

   Int_t n = 0;
   while (pch != NULL)
   {
      fUp[n] = atof(pch);
      n++;
      pch = strtok (NULL, ", ");
   }

   return;
}

//______________________________________________________________________________
void TMFitCollection::SetFitPolMeanE(Char_t* s)
{
   Char_t* pch;
   pch = strtok (s,",");

   Int_t n = 0;
   while (pch != NULL)
   {
      szPolMeanE[n] = pch;
      n++;
      pch = strtok (NULL, ", ");
   }

   return;
}

//______________________________________________________________________________
void TMFitCollection::SetFitPolSigmaE(Char_t* s)
{
   Char_t* pch;
   pch = strtok (s,",");

   Int_t n = 0;
   while (pch != NULL)
   {
      szPolSigmaE[n] = pch;
      n++;
      pch = strtok (NULL, ", ");
   }

   return;
}

//______________________________________________________________________________
void TMFitCollection::SetFitPolLoE(Char_t* s)
{
   Char_t* pch;
   pch = strtok (s,",");

   Int_t n = 0;
   while (pch != NULL)
   {
      szPolLoE[n] = pch;
      n++;
      pch = strtok (NULL, ", ");
   }

   return;
}

//______________________________________________________________________________
void TMFitCollection::SetFitPolUpE(Char_t* s)
{
   Char_t* pch;
   pch = strtok (s,",");

   Int_t n = 0;
   while (pch != NULL)
   {
      szPolUpE[n] = pch;
      n++;
      pch = strtok (NULL, ", ");
   }

   return;
}

//______________________________________________________________________________
TGraphErrors* TMFitCollection::CleanGraphE(Int_t c, TGraphErrors* g)
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
         if (oFit[c][GetEnergyBin(x)]->GetParameter(0)==0.0)
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

   Info("TMFitCollection::CleanGraphE","Removed %i/%i points of %s %s",n,k,g->Class_Name(),g->GetName());

   TGraphErrors* gNew = new TGraphErrors(g->GetN()+2);
   Double_t x = 0.;
   Double_t y = 0.;
   g->GetPoint(0,x,y);
   gNew->SetPoint(0,fMin,y);
   gNew->SetPointError(0,0.,0.);
   g->GetPoint(g->GetN()-1,x,y);
   gNew->SetPoint(gNew->GetN()-1,fMax,y);
   gNew->SetPointError(gNew->GetN()-1,0.,0.);

   for (Int_t i = 0; i < g->GetN(); i++)
   {
      Double_t xn = 0.;
      Double_t yn = 0.;
      g->GetPoint(i,xn,yn);
      Double_t xe = g->GetErrorX(i);
      Double_t ye = g->GetErrorY(i);
      gNew->SetPoint(i+1,xn,yn);
      gNew->SetPointError(i+1,xe,ye);
   }

   return gNew;
}

//______________________________________________________________________________
void TMFitCollection::CleanGraphCT(Int_t e, TGraphErrors* g)
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
         if (oFit[GetCosThetaBin(x)][e]->GetParameter(0)==0.0)
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

   Info("TMFitCollection::CleanGraphCT","Removed %i/%i points of %s %s",n,k,g->Class_Name(),g->GetName());

   return;
}

//______________________________________________________________________________
TH2** TMFitCollection::RecalculateExcFuncErrors(TOEnergyThetaData* eIn, const Char_t* s, Double_t fM, Double_t fB, Double_t fT, Bool_t kA)
{
   TH2** h = new TH2*[fNCTbins];

   for (Int_t i = 0; i < fNCTbins; i++)
   {
      h[i] = (TH2F*)GetExcitationFunction(i);

      Double_t* lo = GetLowerCutGraphE(i)->GetY();
      Double_t* up = GetUpperCutGraphE(i)->GetY();

      for (Int_t j = 0; j < fNEbins; j++)
      {
         TH1D* hSig = (TH1D*)GetSignalObject(i,j)->Clone(Form("hSig_%i_%i",i,j));
         Double_t val = 0.;
         Double_t err_val = 0.;
         Double_t e = eIn->GetValue(eIn->FindEnergyBin(GetEnergy(j)),i)*fB*fT;
//         Double_t err_e = eff->GetError(eff->FindEnergyBin(fEnergy[j]),i)*fBR*fTD;
         if (kAll)
            val = hSig->IntegralAndError(1,hSig->GetNbinsX(),err_val)/e;
         else
            val = hSig->IntegralAndError(hSig->FindBin(lo[j]),hSig->FindBin(up[j]),err_val)/e;
         err_val /= e;//TMath::Sqrt(err_val*err_val/e/e+val*val*err_e*err_e/e/e/e/e);
         if (isnan(val)) val = 0.;
         if (isnan(err_val)) err_val = 0.;
         h[i]->SetBinError(j+1,h[i]->GetYaxis()->FindBin(fM),err_val);
      }
   }

   if (kA)
      Info("TMFitCollection::RecalculateExcFuncErrors","Recalculated Excitation Function Errors over the full range");
   else
      Info("TMFitCollection::RecalculateExcFuncErrors","Recalculated Excitation Function Errors within the cut range");

   return h;
}
//______________________________________________________________________________
TH2** TMFitCollection::RecalculateExcitationFunctions(TOEnergyThetaData* eIn, const Char_t* s, Double_t fM, Double_t fB, Double_t fT, Bool_t kA)
{
   TH2** h = new TH2*[fNCTbins];

   for (Int_t i = 0; i < fNCTbins; i++)
   {
      h[i] = (TH2F*)oSigIn[i]->Clone(Form("%s_%i",s,i));
      for (Int_t j = 0; j < h[i]->GetNbinsX(); j++)
      {
         for (Int_t k = 0; k < h[i]->GetNbinsY(); k++)
         {
            h[i]->SetBinContent(j+1,k+1,0.);
            h[i]->SetBinError(j+1,k+1,0.);
         }
      } 
//      hIM[i] = new TH2F(Form("%s_%i",sExFunc,i),Form("%s_%i",sExFunc,i),
//                  fNEbins,fEnergy[0],fEnergy[fNEbins-1],
//                  hIMref->GetNbinsY(),hIMref->GetYaxis()->GetBinLowEdge(1),hIMref->GetYaxis()->GetBinUpEdge(hIMref->GetNbinsY()));
      h[i]->Sumw2();

      Double_t* lo = GetLowerCutGraphE(i)->GetY();
      Double_t* up = GetUpperCutGraphE(i)->GetY();

      for (Int_t j = 0; j < fNEbins; j++)
      {
         TH1D* hSig = (TH1D*)GetSignalObject(i,j)->Clone(Form("hSig_%i_%i",i,j));
         Double_t val = 0.;
         Double_t err_val = 0.;
         Double_t e = eIn->GetValue(eIn->FindEnergyBin(GetEnergy(j)),i)*fB*fT;
//         Double_t err_e = eff->GetError(eff->FindEnergyBin(fEnergy[j]),i)*fBR*fTD;
         if (kAll)
            val = hSig->IntegralAndError(1,hSig->GetNbinsX(),err_val)/e;
         else
            val = hSig->IntegralAndError(hSig->FindBin(lo[j]),hSig->FindBin(up[j]),err_val)/e;
         err_val /= e;//TMath::Sqrt(err_val*err_val/e/e+val*val*err_e*err_e/e/e/e/e);
         if (isnan(val)) val = 0.;
         if (isnan(err_val)) err_val = 0.;
         h[i]->SetBinContent(j+1,h[i]->GetYaxis()->FindBin(fM),val);
         h[i]->SetBinError(j+1,h[i]->GetYaxis()->FindBin(fM),err_val);
      }
   }

   if (kA)
      Info("TMFitCollection::RecalculateExcitationFunctions","Recalculated Excitation Functions over the full range");
   else
      Info("TMFitCollection::RecalculateExcitationFunctions","Recalculated Excitation Functions within the cut range");

   return h;
}

//______________________________________________________________________________
void TMFitCollection::CalculateExcitationFunctions(Bool_t kEff)
{
   sb = new TOEnergyThetaData(fNEbins, fNCTbins, "sb_ratio");

   for (Int_t i = 0; i < fNCTbins; i++)
   {
      hIM[i] = (TH2F*)oSigIn[i]->Clone(Form("%s_%i",sExFunc,i));
      for (Int_t j = 0; j < hIM[i]->GetNbinsX(); j++)
      {
         for (Int_t k = 0; k < hIM[i]->GetNbinsY(); k++)
         {
            hIM[i]->SetBinContent(j+1,k+1,0.);
            hIM[i]->SetBinError(j+1,k+1,0.);
         }
      } 
//      hIM[i] = new TH2F(Form("%s_%i",sExFunc,i),Form("%s_%i",sExFunc,i),
//                  fNEbins,fEnergy[0],fEnergy[fNEbins-1],
//                  hIMref->GetNbinsY(),hIMref->GetYaxis()->GetBinLowEdge(1),hIMref->GetYaxis()->GetBinUpEdge(hIMref->GetNbinsY()));
      hIM[i]->Sumw2();

//      Double_t* lo = gCutLoE[i]->GetY();
//      Double_t* up = gCutUpE[i]->GetY();

      for (Int_t j = 0; j < fNEbins; j++)
      {
         TH1D* hDat = (TH1D*)oDat[i][j]->Clone(Form("hDat_%i_%i",i,j));
         TH1D* hSig = (TH1D*)oSig[i][j]->Clone(Form("hSig_%i_%i",i,j));
Printf("dat:%f sig:%f",hDat->GetMaximum(),hSig->GetMaximum());
         Double_t val = 0.;
         Double_t vald = 0.;
         Double_t err_val = 0.;
         Double_t err_vald = 0.;
         Double_t e = 1.;
         if (kEff) e = eff->GetValue(eff->FindEnergyBin(fEnergy[j]),i)*fBR*fTD;
//         Double_t err_e = eff->GetError(eff->FindEnergyBin(fEnergy[j]),i)*fBR*fTD;
         Int_t d_lo = 0;
         Int_t d_up = 0;
         Int_t s_lo = 0;
         Int_t s_up = 0;

         if (kAll)
         {
            d_lo = 1;
            d_up = hDat->GetNbinsX();
            s_lo = 1;
            s_up = hSig->GetNbinsX();
         }
         else
         {
            d_lo = hDat->FindBin(fCutLoE[i]->Eval(hIM[i]->GetXaxis()->GetBinCenter(j+1)));
            d_up = hDat->FindBin(fCutUpE[i]->Eval(hIM[i]->GetXaxis()->GetBinCenter(j+1)));
            s_lo = hSig->FindBin(fCutLoE[i]->Eval(hIM[i]->GetXaxis()->GetBinCenter(j+1)));
            s_up = hSig->FindBin(fCutUpE[i]->Eval(hIM[i]->GetXaxis()->GetBinCenter(j+1)));
         }

         vald = hDat->IntegralAndError(d_lo,d_up,err_vald);
         val  = hSig->IntegralAndError(s_lo,s_up,err_val);

         sb->SetEnergy(j, hIM[i]->GetXaxis()->GetBinCenter(j+1));

         if (vald!=0.)
         {
            sb->SetValue(j, i, val/vald);
            sb->SetError(j, i, TMath::Sqrt(err_val*err_val/vald/vald+val*val*err_vald*err_vald/vald/vald/vald/vald));
         }
         else
         {
            sb->SetValue(j, i, 0.);
            sb->SetError(j, i, 0.);
         }
//            val = hSig->IntegralAndError(hSig->FindBin(lo[j]),hSig->FindBin(up[j]),err_val)/e;
Printf("S: %f, S+BG: %f SB: %f", val, vald, val/vald);
         val /= e;
         err_val /= e;//TMath::Sqrt(err_val*err_val/e/e+val*val*err_e*err_e/e/e/e/e);
         if (isnan(val)) val = 0.;
         if (isnan(err_val)) err_val = 0.;
         hIM[i]->SetBinContent(j+1,hIM[i]->GetYaxis()->FindBin(fIM),val);
         hIM[i]->SetBinError(j+1,hIM[i]->GetYaxis()->FindBin(fIM),err_val);
Printf("Counts bin %i (E=%f MeV) Sig:[%f,%f] Dat:[%f,%f]: %f",j+1,hIM[i]->GetXaxis()->GetBinCenter(j+1),hSig->GetBinCenter(s_lo),hSig->GetBinCenter(s_up),hDat->GetBinCenter(d_lo),hDat->GetBinCenter(d_up),val);
Printf("Bins Data:%i Bins MC:%i",hDat->GetNbinsX(),hSig->GetNbinsX());
      }
   }

   if (kAll)
      Info("TMFitCollection::CalculateExcitationFunctions","Calculated Excitation Functions over the full range");
   else
      Info("TMFitCollection::CalculateExcitationFunctions","Calculated Excitation Functions within the cut range");

   return;
}

//______________________________________________________________________________
Int_t TMFitCollection::GetEnergyBin(Double_t f)
{
   if (f < fEnergy[0]){ Printf("too low");
return 0;}
   else if (f >= fEnergy[fNEbins-1]){ Printf("too high");
return fNEbins-1;}
   else
   {
      Int_t n = -1;
      for (Int_t i = 0; i < fNEbins; i++)
      {
         Double_t lo = 0.;
         if (i==0) lo = fEnergy[i];
         else lo = fEnergy[i] - 0.5*(fEnergy[i]-fEnergy[i-1]);
         Double_t up = 0;
         if (i<fNEbins-1) up = fEnergy[i] + 0.5*(fEnergy[i+1]-fEnergy[i]);
         else up = fEnergy[i];

         if (f >= lo && f < up )
         {
            n = i;
            break;
         }
      }
      if (n >= 0) 
         return n;
      else
      {
         Error("TMFitCollection::GetEnergyBin","Energy bin for value %f not found",f);
         return 0;
      }
   }
}

//______________________________________________________________________________
Int_t TMFitCollection::GetCosThetaBin(Double_t f)
{
   if (f < fCTlo[0]) return 0;
   else if (f >= fCTup[fNCTbins-1]) return fNCTbins-1;
   else
   {
      Int_t n = -1;
      for (Int_t i = 0; i < fNCTbins; i++)
      {
         if( f >= fCTlo[i] && f < fCTup[i])
            n = i;
      }
      if (n >= 0) 
         return n;
      else
      {
         Error("TMFitCollection::GetCosThetaBin","cos(theta) bin for value %f not found",f);
         return 0;
      }
   }
}

//______________________________________________________________________________
void TMFitCollection::InitObjects()
{
   oDatIn = new TH2F*[fNCTbins];
   oSigIn = new TH2F*[fNCTbins];
   oBGIn  = new TH2F**[fNCTbins];

   oDat = new TH1D**[fNCTbins];
   oSig = new TH1D**[fNCTbins];
   oBG  = new TH1D***[fNCTbins];
   oBGTot = new TH1D**[fNCTbins];
   oTot = new TH1D**[fNCTbins];
   oFit = new TF1**[fNCTbins];

   fPar    = new Double_t**[fNCTbins];
   fParErr = new Double_t**[fNCTbins];

   for (Int_t i = 0; i < fNCTbins; i++)
   {
      oBGIn[i]  = new TH2F*[fNBG];
      oDat[i]   = new TH1D*[fNEbins];
      oSig[i]   = new TH1D*[fNEbins];
      oBG[i]    = new TH1D**[fNEbins];
      oBGTot[i] = new TH1D*[fNEbins];
      oTot[i]   = new TH1D*[fNEbins];
      oFit[i]   = new TF1*[fNEbins];

      fPar[i]    = new Double_t*[fNEbins];
      fParErr[i] = new Double_t*[fNEbins];
 
      for (Int_t j = 0; j < fNEbins; j++)
      {
         oBG[i][j] = new TH1D*[fNBG];

         fPar[i][j]    = new Double_t[fNBG+1];
         fParErr[i][j] = new Double_t[fNBG+1];

         for (Int_t k = 0; k < fNBG+1; k++)
         {
            fPar[i][j][k]    = 1.;
            fParErr[i][j][k] = 0.;
         }
         oFit[i][j] = new TF1(Form("fGauss_%i_%i",i,j), "gaus", fGmin, fGmax);
      }
   }

   fEnergy = new Double_t[fNEbins];
   for (Int_t j = 0; j < fNEbins; j++)
      fEnergy[j] = 0.;

   fCTlo = new Double_t[fNCTbins];
   fCTup = new Double_t[fNCTbins];
   for (Int_t i = 0; i < fNCTbins; i++)
   {
      fCTlo[i] = 0.;
      fCTup[i] = 0.;
   }

   fFitCT = new TF1**[fNEbins];
   gParCT = new TGraphErrors**[fNEbins];
   for (Int_t i = 0; i < fNEbins; i++)
   {
      fFitCT[i] = new TF1*[fNBG+1];
      gParCT[i] = new TGraphErrors*[fNBG+1];
   }

   fFitE = new TF1**[fNCTbins];
   gParE = new TGraphErrors**[fNCTbins];
   for (Int_t i = 0; i < fNCTbins; i++)
   {
      fFitE[i] = new TF1*[fNBG+1];
      gParE[i] = new TGraphErrors*[fNBG+1];
   }

   fLo = new Double_t[fNCTbins];
   fUp = new Double_t[fNCTbins];
   for (Int_t i = 0; i < fNCTbins; i++)
   {
      fLo[i] = 0.;
      fUp[i] = 0.;
   }

   fCutMeanE  = new TF1*[fNCTbins];
   fCutSigmaE = new TF1*[fNCTbins];
   fCutLoE    = new TF1*[fNCTbins];
   fCutUpE    = new TF1*[fNCTbins];

   gCutMeanE  = new TGraphErrors*[fNCTbins];
   gCutSigmaE = new TGraphErrors*[fNCTbins];
   gCutLoE    = new TGraphErrors*[fNCTbins];
   gCutUpE    = new TGraphErrors*[fNCTbins];
   for (Int_t i = 0; i < fNCTbins; i++)
   {
      gCutMeanE[i]  = new TGraphErrors(fNEbins);
      gCutSigmaE[i] = new TGraphErrors(fNEbins);
      gCutLoE[i]    = new TGraphErrors(fNEbins);
      gCutUpE[i]    = new TGraphErrors(fNEbins);
   }

   fCutMeanCT  = new TF1*[fNEbins];
   fCutSigmaCT = new TF1*[fNEbins];
   fCutLoCT    = new TF1*[fNEbins];
   fCutUpCT    = new TF1*[fNEbins];

   gCutMeanCT  = new TGraphErrors*[fNEbins];
   gCutSigmaCT = new TGraphErrors*[fNEbins];
   gCutLoCT    = new TGraphErrors*[fNEbins];
   gCutUpCT    = new TGraphErrors*[fNEbins];
   for (Int_t i = 0; i < fNEbins; i++)
   {
      gCutMeanCT[i]  = new TGraphErrors(fNCTbins);
      gCutSigmaCT[i] = new TGraphErrors(fNCTbins);
      gCutLoCT[i]    = new TGraphErrors(fNCTbins);
      gCutUpCT[i]    = new TGraphErrors(fNCTbins);
   }

   fChiSquare = new Double_t*[fNCTbins];
   for (Int_t i = 0; i < fNCTbins; i++)
   {
      fChiSquare[i] = new Double_t[fNEbins];
      for (Int_t j = 0; j < fNEbins; j++)
      {
         fChiSquare[i][j] = 999999999.;
      }
   }

   gChiSquareE = new TGraphErrors*[fNCTbins];
   for (Int_t i = 0; i < fNCTbins; i++)
   {
      gChiSquareE[i] = new TGraphErrors(fNEbins);
   }

   gChiSquareCT = new TGraphErrors*[fNEbins];
   for (Int_t i = 0; i < fNEbins; i++)
   {
      gChiSquareCT[i] = new TGraphErrors(fNCTbins);
   }

   szPolMeanE   = new Char_t*[fNCTbins];
   szPolSigmaE  = new Char_t*[fNCTbins];
   szPolLoE     = new Char_t*[fNCTbins];
   szPolUpE     = new Char_t*[fNCTbins];

   eCut = new TOEnergyThetaCut(fNCTbins,"cut");

   hIM = new TH2F*[fNCTbins];

   Info("TMFitCollection::InitObjects","Initialized all objects");

   return;
}

//______________________________________________________________________________
void TMFitCollection::Streamer(TBuffer& R__b)
{
    // Stream a TMFitCollection object.

    UInt_t R__s, R__c;

    if (R__b.IsReading())
    {
        Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
        if (R__v) { }

        // call streamer of parent class
        TNamed::Streamer(R__b);

        // copy data members
        R__b >> fNBG;
        R__b >> fNCTbins;
        R__b >> fNEbins;
        R__b >> fRep;
        R__b >> eCut;
        R__b >> fMin;
        R__b >> fMax;
        R__b >> fGMmin;
        R__b >> fGMmax;
        R__b >> fGMtol;
        R__b >> fGHmin;
        R__b >> fGHmax;
        R__b >> fGSmin;
        R__b >> fGSmax;
        R__b >> sb;

        // copy the objects
        oDatIn  = new TH2F*[fNCTbins];
        oSigIn  = new TH2F*[fNCTbins];
        for (Int_t i = 0; i < fNCTbins; i++)
        {
           oDatIn[i]  = (TH2F*) R__b.ReadObject(TH2F::Class());
           oSigIn[i]  = (TH2F*) R__b.ReadObject(TH2F::Class());
        }

        oBGIn  = new TH2F**[fNCTbins];
        for (Int_t i = 0; i < fNCTbins; i++)
        {
           oBGIn[i] = new TH2F*[fNBG];
           for (Int_t j = 0; j < fNBG; j++)
              oBGIn[i][j] = (TH2F*) R__b.ReadObject(TH2F::Class());
        }

        oDat    = new TH1D**[fNCTbins];
        oSig    = new TH1D**[fNCTbins];
        oBGTot  = new TH1D**[fNCTbins];
        oTot    = new TH1D**[fNCTbins];
        oFit    = new TF1**[fNCTbins];
        for (Int_t i = 0; i < fNCTbins; i++)
        {
           oDat[i]    = new TH1D*[fNEbins];
           oSig[i]    = new TH1D*[fNEbins];
           oBGTot[i]  = new TH1D*[fNEbins];
           oTot[i]    = new TH1D*[fNEbins];
           oFit[i]    = new TF1*[fNEbins];
           for (Int_t j = 0; j < fNEbins; j++)
           {
              oDat[i][j]   = (TH1D*) R__b.ReadObject(TH1D::Class());
              oSig[i][j]   = (TH1D*) R__b.ReadObject(TH1D::Class());
              oBGTot[i][j] = (TH1D*) R__b.ReadObject(TH1D::Class());
              oTot[i][j]   = (TH1D*) R__b.ReadObject(TH1D::Class());
              oFit[i][j]   = (TF1*) R__b.ReadObject(TF1::Class());
           }
        }

        oBG     = new TH1D***[fNCTbins];
        for (Int_t i = 0; i < fNCTbins; i++)
        {
           oBG[i]     = new TH1D**[fNEbins];
           for (Int_t j = 0; j < fNEbins; j++)
           {
              oBG[i][j]    = new TH1D*[fNBG];
              for (Int_t k = 0; k < fNBG; k++)
                 oBG[i][j][k] = (TH1D*) R__b.ReadObject(TH1D::Class());
           }
        }

        fEnergy = new Double_t[fNEbins]; 
        for (Int_t i = 0; i < fNEbins; i++)
           R__b >> fEnergy[i];
 
        fCTlo = new Double_t[fNCTbins];
        fCTup = new Double_t[fNCTbins];
        for (Int_t i = 0; i < fNCTbins; i++)
        {
           R__b >> fCTlo[i];
           R__b >> fCTup[i];
        }

        fFitCT = new TF1**[fNEbins];
        gParCT = new TGraphErrors**[fNEbins];
        for (Int_t i = 0; i < fNEbins; i++)
        {
           fFitCT[i] = new TF1*[fNBG+1];
           gParCT[i] = new TGraphErrors*[fNBG+1];
           for (Int_t j = 0; j < fNBG+1; j++)
           {
              fFitCT[i][j] = (TF1*) R__b.ReadObject(TF1::Class());
              gParCT[i][j] = (TGraphErrors*) R__b.ReadObject(TGraphErrors::Class());
           }
        }
    
        fFitE = new TF1**[fNCTbins];
        gParE = new TGraphErrors**[fNCTbins];
        for (Int_t i = 0; i < fNCTbins; i++)
        {
           fFitE[i] = new TF1*[fNBG+1];
           gParE[i] = new TGraphErrors*[fNBG+1];
           for (Int_t j = 0; j < fNBG+1; j++)
           {
              fFitE[i][j] = (TF1*) R__b.ReadObject(TF1::Class());
              gParE[i][j] = (TGraphErrors*) R__b.ReadObject(TGraphErrors::Class());
           }
        }

        fPar = new Double_t**[fNCTbins];
        fParErr = new Double_t**[fNCTbins];
        for (Int_t i = 0; i < fNCTbins; i++)
        {
           fPar[i] = new Double_t*[fNEbins];
           fParErr[i] = new Double_t*[fNEbins];
           for (Int_t j = 0; j < fNEbins; j++)
           {
              fPar[i][j] = new Double_t[fNBG+1];
              fParErr[i][j] = new Double_t[fNBG+1];
              for (Int_t k = 0; k < fNBG+1; k++)
              {
                 R__b >> fPar[i][j][k];
                 R__b >> fParErr[i][j][k];
              }
           }
        }

        fCutMeanE  = new TF1*[fNCTbins];
        fCutSigmaE = new TF1*[fNCTbins];
        fCutLoE    = new TF1*[fNCTbins];
        fCutUpE    = new TF1*[fNCTbins];
     
        gCutMeanE  = new TGraphErrors*[fNCTbins];
        gCutSigmaE = new TGraphErrors*[fNCTbins];
        gCutLoE    = new TGraphErrors*[fNCTbins];
        gCutUpE    = new TGraphErrors*[fNCTbins];
        for (Int_t i = 0; i < fNCTbins; i++)
        {
           fCutMeanE[i]  = (TF1*) R__b.ReadObject(TF1::Class());
           fCutSigmaE[i] = (TF1*) R__b.ReadObject(TF1::Class());
           fCutLoE[i]    = (TF1*) R__b.ReadObject(TF1::Class());
           fCutUpE[i]    = (TF1*) R__b.ReadObject(TF1::Class());

           gCutMeanE[i]  = (TGraphErrors*) R__b.ReadObject(TGraphErrors::Class());
           gCutSigmaE[i] = (TGraphErrors*) R__b.ReadObject(TGraphErrors::Class());
           gCutLoE[i]    = (TGraphErrors*) R__b.ReadObject(TGraphErrors::Class());
           gCutUpE[i]    = (TGraphErrors*) R__b.ReadObject(TGraphErrors::Class());
        }
     
        fCutMeanCT  = new TF1*[fNEbins];
        fCutSigmaCT = new TF1*[fNEbins];
        fCutLoCT    = new TF1*[fNEbins];
        fCutUpCT    = new TF1*[fNEbins];
     
        gCutMeanCT  = new TGraphErrors*[fNEbins];
        gCutSigmaCT = new TGraphErrors*[fNEbins];
        gCutLoCT    = new TGraphErrors*[fNEbins];
        gCutUpCT    = new TGraphErrors*[fNEbins];
        for (Int_t i = 0; i < fNEbins; i++)
        {
           fCutMeanCT[i]  = (TF1*) R__b.ReadObject(TF1::Class());
           fCutSigmaCT[i] = (TF1*) R__b.ReadObject(TF1::Class());
           fCutLoCT[i]    = (TF1*) R__b.ReadObject(TF1::Class());
           fCutUpCT[i]    = (TF1*) R__b.ReadObject(TF1::Class());

           gCutMeanCT[i]  = (TGraphErrors*) R__b.ReadObject(TGraphErrors::Class());
           gCutSigmaCT[i] = (TGraphErrors*) R__b.ReadObject(TGraphErrors::Class());
           gCutLoCT[i]    = (TGraphErrors*) R__b.ReadObject(TGraphErrors::Class());
           gCutUpCT[i]    = (TGraphErrors*) R__b.ReadObject(TGraphErrors::Class());
        }

        fChiSquare = new Double_t*[fNCTbins];
        for (Int_t i = 0; i < fNCTbins; i++)
        {
           fChiSquare[i] = new Double_t[fNEbins];
           for (Int_t j = 0; j < fNEbins; j++)
           {
              R__b >> fChiSquare[i][j];
           }
        }

        gChiSquareE = new TGraphErrors*[fNCTbins];
        for (Int_t i = 0; i < fNCTbins; i++)
        {
           gChiSquareE[i] = (TGraphErrors*) R__b.ReadObject(TGraphErrors::Class());
        }
     
        gChiSquareCT = new TGraphErrors*[fNEbins];
        for (Int_t i = 0; i < fNEbins; i++)
        {
           gChiSquareCT[i] = (TGraphErrors*) R__b.ReadObject(TGraphErrors::Class());
        }

        hIM     = new TH2F*[fNCTbins];
        for (Int_t i = 0; i < fNCTbins; i++)
        {
           hIM[i] = (TH2F*) R__b.ReadObject(TH2F::Class());
        }

        // check byte count
        R__b.CheckByteCount(R__s, R__c, TMFitCollection::IsA());
    }
    else
    {
        R__c = R__b.WriteVersion(TMFitCollection::IsA(), kTRUE);

        // call streamer of parent class
        TNamed::Streamer(R__b);

        // copy data members
        R__b << fNBG;
        R__b << fNCTbins;
        R__b << fNEbins;
        R__b << fRep;
        R__b << eCut;
        R__b << fMin;
        R__b << fMax;
        R__b << fGMmin;
        R__b << fGMmax;
        R__b << fGMtol;
        R__b << fGHmin;
        R__b << fGHmax;
        R__b << fGSmin;
        R__b << fGSmax;
        R__b << sb;

        // copy the objects

        for (Int_t i = 0; i < fNCTbins; i++)
        {
           R__b.WriteObject(oDatIn[i]);
           R__b.WriteObject(oSigIn[i]);
        }

        for (Int_t i = 0; i < fNCTbins; i++)
           for (Int_t j = 0; j < fNBG; j++)
              R__b.WriteObject(oBGIn[i][j]);

        for (Int_t i = 0; i < fNCTbins; i++)
        {
           for (Int_t j = 0; j < fNEbins; j++)
           {
              R__b.WriteObject(oDat[i][j]);
              R__b.WriteObject(oSig[i][j]);
              R__b.WriteObject(oBGTot[i][j]);
              R__b.WriteObject(oTot[i][j]);
              R__b.WriteObject(oFit[i][j]);
           }
        }

        for (Int_t i = 0; i < fNCTbins; i++)
           for (Int_t j = 0; j < fNEbins; j++)
              for (Int_t k = 0; k < fNBG; k++)
                 R__b.WriteObject(oBG[i][j][k]);

        for (Int_t i = 0; i < fNEbins; i++)
           R__b << fEnergy[i];

        for (Int_t i = 0; i < fNCTbins; i++)
        {
           R__b << fCTlo[i];
           R__b << fCTup[i];
        }

        for (Int_t i = 0; i < fNEbins; i++)
        {
           for (Int_t j = 0; j < fNBG+1; j++)
           {
              R__b.WriteObject(fFitCT[i][j]);
              R__b.WriteObject(gParCT[i][j]);
           }
        }

        for (Int_t i = 0; i < fNCTbins; i++)
        {
           for (Int_t j = 0; j < fNBG+1; j++)
           {
              R__b.WriteObject(fFitE[i][j]);
              R__b.WriteObject(gParE[i][j]);
           }
        }

        for (Int_t i = 0; i < fNCTbins; i++)
        {
           for (Int_t j = 0; j < fNEbins; j++)
           {
              for (Int_t k = 0; k < fNBG+1; k++)
              {
                 R__b << fPar[i][j][k];
                 R__b << fParErr[i][j][k];
              }
           }
        }

        for (Int_t i = 0; i < fNCTbins; i++)
        {
           R__b.WriteObject(fCutMeanE[i]);
           R__b.WriteObject(fCutSigmaE[i]);
           R__b.WriteObject(fCutLoE[i]);
           R__b.WriteObject(fCutUpE[i]);

           R__b.WriteObject(gCutMeanE[i]);
           R__b.WriteObject(gCutSigmaE[i]);
           R__b.WriteObject(gCutLoE[i]);
           R__b.WriteObject(gCutUpE[i]);
        }
     
        for (Int_t i = 0; i < fNEbins; i++)
        {
           R__b.WriteObject(fCutMeanCT[i]);
           R__b.WriteObject(fCutSigmaCT[i]);
           R__b.WriteObject(fCutLoCT[i]);
           R__b.WriteObject(fCutUpCT[i]);

           R__b.WriteObject(gCutMeanCT[i]);
           R__b.WriteObject(gCutSigmaCT[i]);
           R__b.WriteObject(gCutLoCT[i]);
           R__b.WriteObject(gCutUpCT[i]);
        }

        for (Int_t i = 0; i < fNCTbins; i++)
        {
           for (Int_t j = 0; j < fNEbins; j++)
           {
              R__b << fChiSquare[i][j];
           }
        }

        for (Int_t i = 0; i < fNCTbins; i++)
        {
           R__b.WriteObject(gChiSquareE[i]);
        }
     
        for (Int_t i = 0; i < fNEbins; i++)
        {
           R__b.WriteObject(gChiSquareCT[i]);
        }

        for (Int_t i = 0; i < fNCTbins; i++)
           R__b.WriteObject(hIM[i]);

        // set byte count
        R__b.SetByteCount(R__c, kTRUE);
    }
}


//______________________________________________________________________________
