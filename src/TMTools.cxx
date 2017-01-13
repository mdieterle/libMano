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


#include "TMTools.h"

//______________________________________________________________________________
TH1F* TMTools::Get1DHisto(TString strFile, TString strOldName, TString strNewName)
{
  // get a histo "strOldName" from a file "strFile"
  // and get a clone of it "strNewName"

  TH1F* h1;

  TFile *ffile = TFile::Open(strFile);
  if( ffile->IsZombie() )
     Error("TMTools::Get1DHisto","Can't open file %s\n",ffile->GetName());

  h1 = (TH1F*) ffile->Get( strOldName )->Clone( strNewName );

  ffile->Close();

  return h1;
}
//_____________________________________________________________________________
TH2F* TMTools::Get2DHisto(TString strFile, TString strOldName, TString strNewName)
{
  // get the 2D histo "strOldName" from a file "strFile"
  // and return a clone of name "strNewName"

  TH2F* h2;

  TFile *ffile = TFile::Open(strFile);
  if( ffile->IsZombie() )
     Error("TMTools::Get1DHisto","Can't open file %s\n",ffile->GetName());
  
  h2 = (TH2F*) ffile->Get( strOldName )->Clone( strNewName );

  ffile->Close();

  return h2;
}

//_____________________________________________________________________________
//TH3F* TMTools::Get3DHisto(TString strFile, TString strOldName, TString strNewName)
//{
//  // get the 2D histo "strOldName" from a file "strFile"
//  // and return a clone of name "strNewName"
//
//  TH3F* h3;
//
//  TFile *ffile = TFile::Open(strFile);
//  if( ffile->IsZombie() )
//     Error("TMTools::Get1DHisto","Can't open file %s\n",ffile->GetName());
//  
//  h3 = (TH3F*) ffile->Get( strOldName )->Clone( strNewName );
//
//  ffile->Close();
//
//  return h3;
//}
//
//______________________________________________________________________________
void TMTools::Get2DCollectionSize(TFile* f, const Char_t* szColl, Int_t& ia, Int_t& ib)
{
   ia = 0;
   ib = 0;

   TIter next(f->GetListOfKeys());
   TKey *key;
   while ((key=(TKey*)next())) 
   {
      if (!strncmp(key->GetClassName(),"TMObjectCollection",strlen("TMObjectCollection")) && 
          !strncmp(key->GetName(),szColl,strlen(szColl)))
      {
         Int_t aa = 0;
         Int_t bb = 0;
         sscanf(key->GetName(),"%*[^_]_%i_%i",&aa,&bb);
         if (aa == 0) ia++;
         if (bb == 0) ia++;
         if (aa>ia) ia = aa;
         if (bb>ib) ib = bb;
      }
   }

   Info("TMTools::Get2DCollectionSize","Found %ix%i TMObjectCollection Objects",ia,ib);

   return;
}

//______________________________________________________________________________
Double_t TMTools::GetFirstNotEmptyBin(TH1* h)
{
   for (Int_t i = 0; i < h->GetNbinsX(); i++)
      if(h->GetBinContent(i+1) > 0)
         return h->GetXaxis()->GetBinLowEdge(i+1);
   return 0.;
}

//______________________________________________________________________________
Double_t TMTools::GetLastNotEmptyBin(TH1* h)
{
   for (Int_t i = h->GetNbinsX()-1; i>-1; i--)
      if(h->GetBinContent(i+1) > 0)
         return h->GetXaxis()->GetBinUpEdge(i+1);
   return 0.;
}

//______________________________________________________________________________
void TMTools::Draw2DCollection(TFile* f, const Char_t* szColl, Int_t ia, Int_t ib, Int_t iNPadX, Int_t iNPadY, Int_t iNDC, Bool_t kFit, Int_t iRebin)
{
   Char_t szName[256];
   Char_t szzName[256];

   Int_t iNPad = iNPadX*iNPadY;
   Int_t nCTbin = ia;
   Int_t nTCbin = ib;

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
               sprintf(szName, "%s_CT_%i_TC_%i_to_%i",     szColl, i+1, j+1, nTCbin);
               sprintf(szzName,"%s_CT_%i_TC_%i_to_%i.ps", szColl, i+1, j+1, nTCbin);

               Info("TMTools::Draw2DCollection","Created Pads for CTbin %i TCbin %i to %i",i+1,j+1,nTCbin);
            }
            else
            {
               sprintf(szName, "%s_CT_%i_TC_%i_to_%i",     szColl, i+1, j+1, j+iNPad);
               sprintf(szzName,"%s_CT_%i_TC_%i_to_%i.ps", szColl, i+1, j+1, j+iNPad);

               Info("TMTools::Draw2DCollection","Created Pads for CTbin %i TCbin %i to %i",i+1,j+1,j+iNPad);
            }

            c = new TCanvas(szName,szName,1000,1000);
         }

         c->cd();

         sprintf(szName,"IM_CT_%i_TC_%i",i+1,j+1);
         p[j] = new TPad(szName, szName, 0.15+(j%iNPadX)*0.8/iNPadX, 0.9125-((j%iNPad)/iNPadX+1)*0.8/iNPadY, 0.15+(j%iNPadX+1)*0.8/iNPadX, 0.9125-((j%iNPad)/iNPadX)*0.8/iNPadY);

         p[j]->SetLeftMargin(0);
         p[j]->SetRightMargin(0);
         p[j]->SetTopMargin(0);
         p[j]->SetBottomMargin(0);

         p[j]->Draw();
         p[j]->cd();

         sprintf(szName,"%s_%i_%i",szColl,i+1,j+1);
         TMObjectCollection* hc = (TMObjectCollection*)f->Get(szName);

         if (kFit)
//            hc->DrawFitCollection();
            hc->DrawPaperCollection();
         else
            hc->DrawArrayCollection(iRebin);

         if (j%iNPad == iNPad-1 || j == nTCbin-1)
         {
            c->cd();

            TGaxis* aX = 0;

            for (Int_t l = 0; l < iNPadX; l++)
            {
               aX = new TGaxis(0.15+l*0.8/iNPadX,
                                       0.9125-((j%iNPad)/iNPadX+1)*0.8/iNPadY,
                                       0.15+(l+1)*0.8/iNPadX,
                                       0.9125-((j%iNPad)/iNPadX+1)*0.8/iNPadY,
                                       -1000.*0.95,
                                       1000.*0.95,
                                       iNDC,"-");
               aX->SetLabelSize(.03);
               aX->SetLabelOffset(-.03);
               aX->SetLabelFont(42);
               aX->Draw();
            }

            TLine aLine;
            aLine.DrawLineNDC(0.95,0.9125-((j%iNPad)/iNPadX+1)*0.8/iNPadY,0.95,0.9125);

            TLatex title;
            title.SetTextAlign(22);
            title.SetTextSize(0.05);
            title.DrawLatex(0.15+0.4,0.5*(0.1125-aX->GetLabelSize())+(iNPadY-((j%iNPad)/iNPadX+1))*0.8/iNPadY,"#font[41]{M_{X}-M_{N} [MeV]}");
            sprintf(szName,"#font[41]{%2.1f < cos(#theta) < %2.1f}",1.-(i+1)*2./nCTbin,1.-i*2./nCTbin);
            title.DrawLatex(0.15+0.4,0.95,szName);
            title.SetTextSize(0.05);
            title.SetTextAngle(90);
            title.DrawLatex(0.5*0.15,0.9125-0.5*((j%iNPad)/iNPadX+1)/iNPadY*0.8,"#font[41]{Counts [a.u.]}");

            c->Write();

//            c->Print(szzName);
         }
      }
   }

   return;
}

//______________________________________________________________________________
void TMTools::Draw2DCollection(TFile* f, const Char_t* szColl, Int_t ia, Int_t ib, Int_t ic, Int_t iNPadX, Int_t iNPadY, Int_t iNDC, Bool_t kFit, Int_t iRebin)
{
   Char_t szName[256];
   Char_t szzName[256];

   Int_t iNPad = iNPadX*iNPadY;
   Int_t nCTbin = ia;
   Int_t nTCbin = ib;
   Int_t nMMbin = ic;

//   for (Int_t i = 0; i < nCTbin; i++)
   for (Int_t i = 4; i < 5; i++)
   {
      for (Int_t j = 0; j < nTCbin; j++)
//      for (Int_t j = 25; j < 26; j++)
      {
         TCanvas* c = 0;
         TPad* p[nMMbin];

         for (Int_t k = 0; k < nMMbin; k++)
         {
            if (k%iNPad == 0)
            {
               if (k+iNPad>nMMbin)
               {
                  sprintf(szName, "%s_CT_%i_TC_%i_MM_%i_to_%i",    szColl, i+1, j+1, k+1, nMMbin);
                  sprintf(szzName,"%s_CT_%i_TC_%i_MM_%i_to_%i.ps", szColl, i+1, j+1, k+1, nMMbin);

                  Info("TMTools::Draw2DCollection","Created Pads for CTbin %i TCbin %i MMbin %i to %i",i+1,j+1,k+1,nMMbin);
               }
               else
               {
                  sprintf(szName, "%s_CT_%i_TC_%i_MM_%i_to_%i",    szColl, i+1, j+1, k+1, k+iNPad);
                  sprintf(szzName,"%s_CT_%i_TC_%i_MM_%i_to_%i.ps", szColl, i+1, j+1, k+1, k+iNPad);

                  Info("TMTools::Draw2DCollection","Created Pads for CTbin %i TCbin %i MMbin %i to %i",i+1,j+1,k+1,k+iNPad);
               }

               c = new TCanvas(szName,szName,1000,1000);
            }

            c->cd();

            sprintf(szName,"IM_CT_%i_TC_%i_MM_%i",i+1,j+1,k+1);
            p[k] = new TPad(szName, szName, 0.15+(k%iNPadX)*0.8/iNPadX, 0.9125-((k%iNPad)/iNPadX+1)*0.8/iNPadY, 0.15+(k%iNPadX+1)*0.8/iNPadX, 0.9125-((k%iNPad)/iNPadX)*0.8/iNPadY);

            p[k]->SetLeftMargin(0);
            p[k]->SetRightMargin(0);
            p[k]->SetTopMargin(0);
            p[k]->SetBottomMargin(0);

            p[k]->Draw();
            p[k]->cd();
            sprintf(szName,"%s_%i_%i_%i",szColl,i+1,j+1,k+1);

            TMObjectCollection* hc = (TMObjectCollection*)f->Get(szName);

            if (kFit)
               hc->DrawFitCollection(0.,300.);
            else
               hc->DrawArrayCollection(iRebin);

            if (k%iNPad == iNPad-1 || k == nMMbin-1)
            {
               c->cd();

               TGaxis* aX = 0;

               for (Int_t l = 0; l < iNPadX; l++)
               {
                  aX = new TGaxis(0.15+l*0.8/iNPadX,
                                          0.9125-((k%iNPad)/iNPadX+1)*0.8/iNPadY,
                                          0.15+(l+1)*0.8/iNPadX,
                                          0.9125-((k%iNPad)/iNPadX+1)*0.8/iNPadY,
                                          0.*0.95,
                                          300.*0.95,
                                          iNDC,"-");
                  aX->SetLabelSize(.03);
                  aX->SetLabelOffset(-.03);
                  aX->SetLabelFont(42);
                  aX->Draw();
               }

               TLine aLine;
               aLine.DrawLineNDC(0.95,0.9125-((k%iNPad)/iNPadX+1)*0.8/iNPadY,0.95,0.9125);

               TLatex title;
               title.SetTextAlign(22);
               title.SetTextSize(0.05);
               title.DrawLatex(0.15+0.4,0.5*(0.1125-aX->GetLabelSize())+(iNPadY-((k%iNPad)/iNPadX+1))*0.8/iNPadY,"#font[41]{M_{#gamma#gamma} [MeV]}");
               sprintf(szName,"#font[41]{%2.1f < cos(#theta) < %2.1f}",1.-(i+1)*2./nCTbin,1.-i*2./nCTbin);
               title.DrawLatex(0.15+0.4,0.95,szName);
               title.SetTextSize(0.05);
               title.SetTextAngle(90);
               title.DrawLatex(0.5*0.15,0.9125-0.5*((k%iNPad)/iNPadX+1)/iNPadY*0.8,"#font[41]{Counts [a.u.]}");

               c->Write();

//               c->Print(szzName);
            }
         }
      }
   }

   return;
}
//______________________________________________________________________________
void TMTools::Draw3DCollection(const Char_t* szFile, const Char_t* szColl, Int_t iNPadX, Int_t iNPadY, Int_t iNDC)
{
   Char_t szName[256];
   Char_t szzName[256];

   Int_t iNPad = iNPadX*iNPadY;

   sprintf(szName,"%s.root",szFile);
   TFile* f = new TFile(szName);
   sprintf(szName,"Plot_%s.root",szFile);
   TFile* g = new TFile(szName,"recreate");

   sprintf(szName,"%s_1_1_1",szColl);
   Int_t nCTbin = ((TMObjectCollection*)f->Get(szName))->GetNCTBins();
   Int_t nTCbin = ((TMObjectCollection*)f->Get(szName))->GetNTCBins();
   Int_t nMMbin = ((TMObjectCollection*)f->Get(szName))->GetNMMBins();

//printf("CT: %i TC: %i MM: %i\n",nCTbin,nTCbin,nMMbin);

   for (Int_t i = 0; i < nCTbin; i++)
   {
      for (Int_t j = 0; j < nTCbin; j++)
      {
         TCanvas* c = 0;
         TPad* p[nMMbin];

         for (Int_t k = 0; k < nMMbin; k++)
         {
//            printf("outside k: %i\n",k);
            if (k%iNPad == 0)
            {
//            printf("inside k: %i\n",k);
               if (k+iNPad>nMMbin)
               {
                  sprintf(szName, "IM_CT_%i_TC_%i_MM_%i_to_%i",     i+1, j+1, k+1, nMMbin);
                  sprintf(szzName,"IM_CT_%i_TC_%i_MM_%i_to_%i.eps", i+1, j+1, k+1, nMMbin);

                  printf("Created Pads for CTbin %i TCbin %i and MMbin %i to %i\n",i+1,j+1,k+1,nMMbin);
               }
               else
               {
                  sprintf(szName, "IM_CT_%i_TC_%i_MM_%i_to_%i",     i+1, j+1, k+1, k+iNPad);
                  sprintf(szzName,"IM_CT_%i_TC_%i_MM_%i_to_%i.eps", i+1, j+1, k+1, k+iNPad);

                  printf("Created Pads for CTbin %i TCbin %i and MMbin %i to %i\n",i+1,j+1,k+1,k+iNPad);
               }

               c = new TCanvas(szName,szName,1000,1000);
            }

            c->cd();

            sprintf(szName,"IM_CT_%i_TC_%i_MM_%i",i+1,j+1,k+1);
            p[k] = new TPad(szName, szName, 0.15+(k%iNPadX)*0.8/iNPadX, 0.9125-((k%iNPad)/iNPadX+1)*0.8/iNPadY, 0.15+(k%iNPadX+1)*0.8/iNPadX, 0.9125-((k%iNPad)/iNPadX)*0.8/iNPadY);

            p[k]->SetLeftMargin(0);
            p[k]->SetRightMargin(0);
            p[k]->SetTopMargin(0);
            p[k]->SetBottomMargin(0);

            p[k]->Draw();
            p[k]->cd();

            sprintf(szName,"%s_%i_%i_%i",szColl,i+1,j+1,k+1);
            TMObjectCollection* hc = (TMObjectCollection*)f->Get(szName);
            hc->Draw();
//            hc->GetBackgroundHisto()->Draw();
//            hc->GetSignal1Histo()->Draw("same");
//            hc->GetSignal2Histo()->Draw("same");
//            hc->GetTotalHisto()->Draw("same");
//            hc->GetDataHisto()->Draw("same");

            if (k%iNPad == iNPad-1 || k == nMMbin-1)
            {
               c->cd();
   
               TGaxis* aX = 0;
   
               for (Int_t l = 0; l < iNPadX; l++)
               {
                  aX = new TGaxis(0.15+l*0.8/iNPadX,
                                          0.9125-((k%iNPad)/iNPadX+1)*0.8/iNPadY,
                                          0.15+(l+1)*0.8/iNPadX,
                                          0.9125-((k%iNPad)/iNPadX+1)*0.8/iNPadY,
                                          0.*0.95,
                                          300.*0.95,
//                                          ((TH1*)hc->GetBackgroundHisto())->GetXaxis()->GetBinLowEdge(1)*0.95,
//                                          ((TH1*)hc->GetBackgroundHisto())->GetXaxis()->GetBinUpEdge(((TH1*)hc->GetBackgroundHisto())->GetNbinsX())*0.95,
                                          iNDC,"-");
                  aX->SetLabelSize(.03);
                  aX->SetLabelOffset(-.03);
                  aX->SetLabelFont(42);
                  aX->Draw();
               }
   
               TLine aLine;
               aLine.DrawLineNDC(0.95,0.9125-((k%iNPad)/iNPadX+1)*0.8/iNPadY,0.95,0.9125);
   
               TLatex title;
               title.SetTextAlign(22);
               title.SetTextSize(0.05);
               title.DrawLatex(0.15+0.4,0.5*(0.1125-aX->GetLabelSize())+(iNPadY-((k%iNPad)/iNPadX+1))*0.8/iNPadY,"#font[41]{M_{#gamma#gamma} [MeV]}");
               title.SetTextSize(0.05);
               title.SetTextAngle(90);
               title.DrawLatex(0.5*0.15,0.9125-0.5*((k%iNPad)/iNPadX+1)/iNPadY*0.8,"#font[41]{Counts [a.u.]}");

               g->cd();
               c->Write();

//               delete c;
//               delete aX;
            }

//            delete [] p;
         }
      }
   }

   f->Close();
   g->Close();

   return;
}

//______________________________________________________________________________
void TMTools::GetSBfromHistoCollection(const Char_t* szFile, const Char_t* szMM, const Char_t* szColl)
{
   Char_t szName[256];
   Char_t szzName[256];

   sprintf(szName,"%s.root",szFile);
   TFile* f = new TFile(szName);
   sprintf(szName,"MM_SB_%s.root",szFile);
   TFile* g = new TFile(szName,"recreate");

   sprintf(szName,"%s_1_1_1",szColl);
   Int_t nCTbin = ((TMObjectCollection*)f->Get(szName))->GetNCTBins();
   Int_t nTCbin = ((TMObjectCollection*)f->Get(szName))->GetNTCBins();
   Int_t nMMbin = ((TMObjectCollection*)f->Get(szName))->GetNMMBins();

   for (Int_t i = 0; i < nCTbin; i++)
   {
      for (Int_t j = 0; j < nTCbin; j++)
      {
         sprintf(szName,"%s_%i_%i",szMM,i+1,j+1);
         sprintf(szzName,"hMM_%i_%i",i+1,j+1);
         TH1D* hMM = (TH1D*)f->Get(szName)->Clone(szzName);
         sprintf(szzName,"hMMsb_%i_%i",i+1,j+1);
         TH1D* hMMsb = (TH1D*)f->Get(szName)->Clone(szzName);
         sprintf(szName,"hSB_%i_%i",i+1,j+1);
         TH1D* hSB = new TH1D(szName,szName,hMM->GetNbinsX(),hMM->GetXaxis()->GetBinLowEdge(1),hMM->GetXaxis()->GetBinUpEdge(hMM->GetNbinsX()));

         for (Int_t k = 0; k < nMMbin; k++)
         {
            sprintf(szName,"%s_%i_%i_%i",szColl,i+1,j+1,k+1);
            TMObjectCollection* hc = (TMObjectCollection*)f->Get(szName);

            Double_t fSige = 0.;
            Double_t fTote = 0.;

            Double_t fSig = ((TH1*)hc->GetSignalObject())->IntegralAndError(1,((TH1*)hc->GetSignalObject())->GetNbinsX(),fSige);
            Double_t fTot = ((TH1*)hc->GetTotalObject())->IntegralAndError(1,((TH1*)hc->GetTotalObject())->GetNbinsX(),fTote);

            Double_t fSB = 0.;
            Double_t fSBe = 0.;

            if (fTot == 0. || fSig == 0.)
            {
               fSB = 0.;
               fSBe = 0.;
            }
            else
            {
               fSB = fSig/fTot;
               fSBe = TMath::Sqrt(TMath::Power(fSig/fTot/fTot*fTote,2.) + TMath::Power(fSige/fTot,2.));
            }

//            printf("CTbin %i TCbin %i MMbin %i - SB %5.3f SBe: %5.3f\n",i+1,j+1,k+1,fSB,fSBe);

            hSB->SetBinContent(k+1,fSB);
            hSB->SetBinError(k+1,fSBe);

            delete hc;
         }

         hMMsb->Multiply(hSB);
         g->cd();
         hMM->Write();
         hSB->Write();
         hMMsb->Write();

         delete hSB;
         delete hMM;
         delete hMMsb;

         printf("Corrected MM histo of CTbin %i and TCbin %i\n",i+1,j+1);
      }
   }

   f->Close();
   g->Close();

   return;
}

//______________________________________________________________________________
void TMTools::GraphToPolygonRight(TFile* f, const Char_t* szGraph, Double_t fLeft, Double_t fBottom, Double_t fRight, Double_t fTop, Double_t fSup)
{
   TGraphErrors* g = (TGraphErrors*)f->Get(szGraph)->Clone("g"); 
   Int_t ng = g->GetN()-2;

   Double_t xg[ng];
   Double_t yg[ng];

   for (Int_t i = 2; i < ng+2; i++)
//   {
      g->GetPoint(i,xg[i-2],yg[i-2]);
//      Info("TMTools::GraphToPolygonRight","Found Point %i[%5.3f,%5.3f]",i+1,xg[i],yg[i]);
//   }

   Int_t na = 6;
   Int_t n = ng+na;

   TCutG* cut = new TCutG("cut",n);

   Double_t x[na];
   Double_t y[na];

   x[0] = xg[ng-1];
   y[0] = yg[ng-1];
   x[1] = xg[ng-1];
   y[1] = fTop;
   x[2] = fLeft;
   y[2] = fTop;
   x[3] = fLeft;
   y[3] = fBottom;
   x[4] = fRight;
   y[4] = fBottom;
   x[5] = fRight;
   if (fSup == 0.)
      y[5] = yg[0];
   else
      y[5] = fSup;

//   Double_t* xn;
//   Double_t* yn;
//
//   xn = &x[0];
//   yn = &y[0];
//
//   for (Int_t i = na; i<n;i++)
//   {
//      xn[i] = xg[i-na];
//      yn[i] = yg[i-na];
//   }

   for (Int_t i = 0; i < n; i++)
   {
      if (i<na)
      {
         cut->SetPoint(i,x[i],y[i]);
         Info("TMTools::GraphToPolygonRight","Added Point %i[%5.3f,%5.3f]",i+1,x[i],y[i]);
      }

//         Info("TMTools::GraphToPolygonRight","Added Point %i[%5.3f,%5.3f] - [%5.3f,%5.3f]",i+1,xn[i],yn[i],x[i],y[i]);
      else
      {
         cut->SetPoint(i,xg[i-na],yg[i-na]);
         Info("TMTools::GraphToPolygonRight","Added Point %i[%5.3f,%5.3f]",i+1,xg[i-na],yg[i-na]);
      }
//         Info("TMTools::GraphToPolygonRight","Added Point %i[%5.3f,%5.3f] - [%5.3f,%5.3f]",i+1,xn[i],yn[i],xg[i-na],yg[i-na]);

   }


//   TCutG* cut = new TCutG("cut",n,xn,yn);

   cut->Write();

}

//______________________________________________________________________________
void TMTools::GraphToPolygonTop(TFile* f, const Char_t* szGraph, Double_t fLeft, Double_t fBottom, Double_t fRight)
{
   TGraphErrors* g = (TGraphErrors*)f->Get(szGraph)->Clone("g"); 
   Int_t ng = 61;//g->GetN()-2;

   Double_t xg[ng];
   Double_t yg[ng];

   for (Int_t i = 2; i < ng+2; i++)
   {
      g->GetPoint(i,xg[i-2],yg[i-2]);
      Info("TMTools::GraphToPolygonRight","Found Point %i[%5.3f,%5.3f]",i,xg[i-2],yg[i-2]);
   }

   Int_t na = 5;
   Int_t n = ng+na;

   TCutG* cut = new TCutG("cut",n);

   Double_t x[na];
   Double_t y[na];

   x[0] = xg[ng-1];
   y[0] = yg[ng-1];
   x[1] = fRight;
   y[1] = yg[ng-1];
   x[2] = fRight;
   y[2] = fBottom;
   x[3] = fLeft;
   y[3] = fBottom;
   x[4] = fLeft;
   y[4] = yg[0];

//   Double_t* xn;
//   Double_t* yn;
//
//   xn = &x[0];
//   yn = &y[0];
//
//   for (Int_t i = na; i<n;i++)
//   {
//      xn[i] = xg[i-na];
//      yn[i] = yg[i-na];
//   }

   for (Int_t i = 0; i < n; i++)
   {
      if (i<na)
      {
         cut->SetPoint(i,x[i],y[i]);
         Info("TMTools::GraphToPolygonRight","Added Point %i[%5.3f,%5.3f]",i+1,x[i],y[i]);
      }

//         Info("TMTools::GraphToPolygonRight","Added Point %i[%5.3f,%5.3f] - [%5.3f,%5.3f]",i+1,xn[i],yn[i],x[i],y[i]);
      else
      {
         cut->SetPoint(i,xg[i-na],yg[i-na]);
         Info("TMTools::GraphToPolygonRight","Added Point %i[%5.3f,%5.3f]",i+1,xg[i-na],yg[i-na]);
      }
//         Info("TMTools::GraphToPolygonRight","Added Point %i[%5.3f,%5.3f] - [%5.3f,%5.3f]",i+1,xn[i],yn[i],xg[i-na],yg[i-na]);

   }


//   TCutG* cut = new TCutG("cut",n,xn,yn);

   cut->Write();

}


//______________________________________________________________________________
void TMTools::GetMeanAndSigma(TFile* f, const Char_t* szFunc, Int_t ia, Int_t ib)
{
  Char_t szName[256];

  for (Int_t i = 0; i < ia; i++)
  {
     TGraphErrors* gMean  = new TGraphErrors(ib);
     sprintf(szName,"gMean_%i",i+1);
     gMean->SetName(szName);
     TGraphErrors* gSigma = new TGraphErrors(ib);
     sprintf(szName,"gSigma_%i",i+1);
     gSigma->SetName(szName);

     for (Int_t j = 0; j < ib; j++)
     {
        sprintf(szName,"%s_%i_%i",szFunc,i+1,j+1);
printf("%s\n",szName);
        TF1* ff = (TF1*)f->Get(szFunc)->Clone("ff");
        gMean->SetPoint(i,atof(ff->GetTitle()),ff->GetParameter(1));
        gMean->SetPointError(i,0.,ff->GetParError(1));
        gSigma->SetPoint(i,atof(ff->GetTitle()),ff->GetParameter(2));
        gSigma->SetPointError(i,0.,ff->GetParError(2));
        delete ff;
     }
     
     gMean->Write(); 
     gSigma->Write();
     delete gMean;
     delete gSigma;
  }

   return;
}

//______________________________________________________________________________
void TMTools::DrawCutQuality(TFile* f, TFile* fData, const Char_t* szHisto, Int_t iCT, 
                             Double_t fSigma, const Char_t* szPol, Double_t fLoLeft, Double_t fLoRight,
                             Double_t fUpLeft, Double_t fUpRight)
{
   Char_t szName[256];
   Char_t szzName[256];

   TOEnergyThetaCut* eC = (TOEnergyThetaCut*)f->Get("mm_cut")->Clone("eC");
   TOEnergyThetaCut* mC = new TOEnergyThetaCut(iCT,"mm_cut_func");  

   Double_t fMax = 1500.;
   if (fLoRight>1500.)
      fMax=2500.;

   if (!fUpLeft) fUpLeft = fLoLeft;
   if (!fUpRight) fUpRight = fLoRight;

   for (Int_t i = 0; i < iCT; i++) 
   {
      sprintf(szName,"%s_%i",szHisto,i);
      sprintf(szzName,"hMM_%i",i);
      TH2F* h = (TH2F*)fData->Get(szName)->Clone(szzName);

      TGraphErrors* gM = eC->GetMeanGraph(i);
      TGraphErrors* gS = eC->GetSigmaGraph(i);
      TGraphErrors* gL = new TGraphErrors(gM->GetN());
      TGraphErrors* gR = new TGraphErrors(gM->GetN());
      gL->SetMarkerStyle(20);
      gL->SetMarkerColor(2);
      gR->SetMarkerStyle(20);
      gR->SetMarkerColor(2);

      for (Int_t j = 0; j < gM->GetN(); j++)
      {
         Double_t xm = 0.;
         Double_t xs = 0.;
         Double_t ym = 0.;
         Double_t ys = 0.;

         gM->GetPoint(j,xm,ym);
         gS->GetPoint(j,xs,ys);

//printf   ("CT: %i Energy: %5.3f  |  Left Value: %5.3f  Mean: %5.3f  Right value: %5.3f\n", i, xm, ym-fSigma*ys, ym, ym+fSigma*ys);

         gL->SetPoint(j,xm,ym-fSigma*ys);
         gR->SetPoint(j,xm,ym+fSigma*ys);
         gL->SetPointError(j,gM->GetErrorX(j),gM->GetErrorY(j)+fSigma*gS->GetErrorY(j));
         gR->SetPointError(j,gM->GetErrorX(j),gM->GetErrorY(j)+fSigma*gS->GetErrorY(j));
      }

      sprintf(szName,"MM_low_%i",i);
      TF1* fL = new TF1(szName,szPol,400,fMax);
      sprintf(szName,"MM_up_%i",i);
      TF1* fR = new TF1(szName,szPol,400,fMax);

      fL->SetLineWidth(2);
      fL->SetLineColor(2);
      fR->SetLineWidth(2);
      fR->SetLineColor(2);

//      if (fLoRight>1500.)
//      {
//         if (i == 0)
//            fLoLeft = 650.;
//         if (i == 1)
//            fLoLeft = 600.;
//         if (i == 2)
//            fLoRight = 1800.;
//         fUpLeft=fLoLeft;
//      }

      gL->Fit(fL,"+0MQ","",fLoLeft,fLoRight);
      gR->Fit(fR,"+0MQ","",fUpLeft,fUpRight);

      Int_t nPol = 0;
      sscanf(szPol,"pol%d",&nPol);
      Int_t nPar = nPol+1;
      Int_t nDim = 2.*nPar;

      Double_t par[nDim];
      fR->GetParameters(&par[0]);
      fL->GetParameters(&par[nPar]);

      sprintf(szName,"MM_mean_%i",i);
      TF1* fM = new TF1(szName,Form("(%s(0)+%s(%d))/2.",szPol,szPol,nPar),400,fMax);
      sprintf(szName,"MM_sigma_%i",i);
      TF1* fS = new TF1(szName,Form("(%s(0)-%s(%d))/2./%lf",szPol,szPol,nPar,fSigma),400,fMax);

      fM->SetParameters(par);
      fS->SetParameters(par);

      sprintf(szName,"MM_Eg_CT_%i",i);
      TCanvas* c = new TCanvas(szName,szName,800,600);
      c->cd();
      h->Draw("col");
      gL->Draw("same"); 
      gR->Draw("same"); 
      fL->Draw("same");
      fR->Draw("same");

      mC->SetLowerGraph(i,gL);
      mC->SetUpperGraph(i,gR);
      mC->SetMeanGraph(i,gM);
      mC->SetSigmaGraph(i,gS);

      mC->SetLowerFunction(i,fL);
      mC->SetUpperFunction(i,fR);
      mC->SetMeanFunction(i,fM);
      mC->SetSigmaFunction(i,fS);

      c->Write();
      delete c;
      delete gL;
      delete gR;
      delete fL;
      delete fR;
      delete fM;
      delete fS;
   }

   mC->Write();
   delete mC;

   return;
}

//______________________________________________________________________________
void TMTools::DrawCopCutQuality(TFile* f, TFile* fData, const Char_t* szHisto, Int_t iCT, 
                             Double_t fSigma, const Char_t* szPol, Double_t fLeft, Double_t fRight)
{
   Char_t szName[256];
   Char_t szzName[256];

   TOEnergyThetaCut* eC = (TOEnergyThetaCut*)f->Get("cop_cut")->Clone("eC");
   TOEnergyThetaCut* mC = new TOEnergyThetaCut(iCT,"cop_cut_func");  

   Double_t fMax = 1500.;
   if (fRight>1500.)
      fMax=2500.;

//   Double_t parl[3];
//   Double_t parr[3];

   for (Int_t i = 0; i < iCT; i++) 
   {
      sprintf(szName,"%s_%i",szHisto,i);
      sprintf(szzName,"hMM_%i",i);
      TH2F* h = (TH2F*)fData->Get(szName)->Clone(szzName);

      TGraphErrors* gM = eC->GetMeanGraph(i);
      TGraphErrors* gS = eC->GetSigmaGraph(i);
      TGraphErrors* gL = new TGraphErrors(gM->GetN());
      TGraphErrors* gR = new TGraphErrors(gM->GetN());
      gL->SetMarkerStyle(20);
      gL->SetMarkerColor(2);
      gR->SetMarkerStyle(20);
      gR->SetMarkerColor(2);

      for (Int_t j = 0; j < gM->GetN(); j++)
      {
         Double_t xm = 0.;
         Double_t xs = 0.;
         Double_t ym = 0.;
         Double_t ys = 0.;

         gM->GetPoint(j,xm,ym);
         gS->GetPoint(j,xs,ys);

printf   ("CT: %i Energy: %5.3f  |  Left Value: %5.3f  Mean: %5.3f  Right value: %5.3f\n", i, xm, ym-fSigma*ys, ym, ym+fSigma*ys);

         gL->SetPoint(j,xm,ym-fSigma*ys);
         gR->SetPoint(j,xm,ym+fSigma*ys);
         gL->SetPointError(j,gM->GetErrorX(j),gM->GetErrorY(j)+fSigma*gS->GetErrorY(j));
         gR->SetPointError(j,gM->GetErrorX(j),gM->GetErrorY(j)+fSigma*gS->GetErrorY(j));
      }

      sprintf(szName,"MM_low_%i",i);
      TF1* fL = new TF1(szName,szPol,400,fMax);
      sprintf(szName,"MM_up_%i",i);
      TF1* fR = new TF1(szName,szPol,400,fMax);

      fL->SetLineWidth(2);
      fL->SetLineColor(2);
      fR->SetLineWidth(2);
      fR->SetLineColor(2);

      if (i == iCT - 1)
      {
         fLeft = 750.;
         if (fRight>1800.)
            fRight = 2300.;
         else
            fRight = 1400.;
      }

//      Int_t nparl = fL->GetNpar();
//      Double_t parl[nparl];
//      Int_t nparr = fR->GetNpar();
//      Double_t parr[nparr];

//      if (i != iCT-1)
//      {
       gL->Fit(fL,"+R0MQ","",fLeft,fRight);
       gR->Fit(fR,"+R0MQ","",fLeft,fRight);

//         fL->GetParameters(&parl[0]);
//         fR->GetParameters(&parr[0]);
//      }
//      else
//      {
//         fL->SetParameters(parl);
//         fR->SetParameters(parr);
//      }

      sprintf(szName,"MM_Eg_CT_%i",i);
      TCanvas* c = new TCanvas(szName,szName,800,600);
      c->cd();
      h->Draw("col");
      gL->Draw("same"); 
      gR->Draw("same"); 
      fL->Draw("same");
      fR->Draw("same");

      mC->SetLowerFunction(i,fL);
      mC->SetUpperFunction(i,fR);

      c->Write();
//      fL->Write();
//      fR->Write();
      delete c;
      delete gL;
      delete gR;
//      delete fL;
//      delete fR;
   }

   mC->Write();
   delete mC;

   return;
}

//______________________________________________________________________________
void TMTools::CreateExcitationFunction(TFile* f, const Char_t* szObj, TFile* fEff, Int_t nCTbin, Int_t nTCbin, const Char_t* szHisto, Double_t nBR, Double_t nTD)
{
   Char_t szName[256];

   TH2F** fH_ExFunc;
   fH_ExFunc = new TH2F*[nCTbin];

   TOEnergyThetaData* eff = (TOEnergyThetaData*)fEff->Get("eff")->Clone("eff");
   TOEnergyThetaCut* mm = (TOEnergyThetaCut*)f->Get("mm_cut_func")->Clone("mm");

   for (Int_t i = 0; i < nCTbin; i++)
   {
      sprintf(szName, "%s_%d", szHisto, i);
//      fH_ExFunc[i] = new TH2F(szName,szName,55,390,1515,400,0,1000);
      fH_ExFunc[i] = new TH2F(szName,szName,40,400,1600,500,0,1000);
      fH_ExFunc[i]->Sumw2();

      TF1* fLo = mm->GetLowerFunction(i);
      TF1* fUp = mm->GetUpperFunction(i);

      for (Int_t j = 0; j < nTCbin; j++)
      {
         Double_t e = eff->GetValue(j, i);
         if (e < 1e-6) e = 1;

         sprintf(szName,"%s_%i_%i",szObj,i+1,j+1);
         TMObjectCollection* hc = (TMObjectCollection*)f->Get(szName);

         sprintf(szName,"hSig_%i_%i",i+1,j+1);
         TH1D* hSig = (TH1D*)hc->GetSignalObject()->Clone(szName);

         Double_t err = 0.;

         Double_t xlo = fLo->Eval(fH_ExFunc[i]->GetXaxis()->GetBinCenter(j+1));
         Double_t xup = fUp->Eval(fH_ExFunc[i]->GetXaxis()->GetBinCenter(j+1));

         fH_ExFunc[i]->SetBinContent(j+1,eff->GetNEnergyBin(),hSig->IntegralAndError(hSig->FindBin(xlo),hSig->FindBin(xup),err) / e / nBR / nTD);
         fH_ExFunc[i]->SetBinError(j+1,eff->GetNEnergyBin(),err / e);
      }

      fH_ExFunc[i]->Write();
   }

   return;
}

//______________________________________________________________________________
void TMTools::CorrectEmptyTarget(TFile* fCS, TFile* fCSempty, Double_t scale)
{
    TOEnergyThetaData* cs  = (TOEnergyThetaData*)fCS->Get("cs")->Clone("cs");
    TOEnergyThetaData* cse = (TOEnergyThetaData*)fCSempty->Get("cs")->Clone("cse");

    if (cs->GetNEnergyBin() != cse->GetNEnergyBin())
       Error("TMTools::CorrectEmptyTarget","Number of energy bins (%i) in file %s differs from number (%i) in file %s",
             cs->GetNEnergyBin(),fCS->GetName(),cse->GetNEnergyBin(),fCSempty->GetName());
    if (cs->GetNCosThetaBin() != cse->GetNCosThetaBin())
       Error("TMTools::CorrectEmptyTarget","Number of CT bins (%i) in file %s differs from number (%i) in file %s",
             cs->GetNCosThetaBin(),fCS->GetName(),cse->GetNCosThetaBin(),fCSempty->GetName());

    cse->GetData()->Scale(scale);

    TOEnergyThetaData* csc = new TOEnergyThetaData(cs->GetNEnergyBin(),cs->GetNCosThetaBin(),"cs");

    for (Int_t i = 0; i < cs->GetNEnergyBin(); i++)
    {
        // set energy
        csc->SetEnergy(i, cs->GetEnergy(i));
        csc->SetEnergyError(i, cs->GetEnergyError(i));

        for (Int_t j = 0; j < cs->GetNCosThetaBin(); j++)
        {
            Double_t v_cs = cs->GetValue(i,j);
            Double_t ve_cs = cs->GetError(i,j);
            Double_t v_cse = cse->GetValue(i,j);
            Double_t ve_cse = cse->GetError(i,j);

            csc->SetValue(i, j, v_cs-v_cse);
            csc->SetError(i, j, TMath::Sqrt(ve_cs*ve_cs+ve_cse*ve_cse));
        }
    }

    csc->Write();
}

//______________________________________________________________________________
TOEnergyThetaData* TMTools::AverageTOEnergyThetaData(TOEnergyThetaData* e1, TOEnergyThetaData* e2)
{
    Int_t iEbin = e1->GetNEnergyBin();
    if (e2->GetNEnergyBin() > iEbin) iEbin = e2->GetNEnergyBin();

    Info("TMTools::AverageTOEnergyThetaData","Creating %i energy bins from data 1 (%i) and data 2 (%i)",
         iEbin,e1->GetNEnergyBin(),e2->GetNEnergyBin());

    Int_t iCTbin = e1->GetNCosThetaBin();
    if (e2->GetNCosThetaBin() > iCTbin) iCTbin = e2->GetNCosThetaBin();

    Info("TMTools::AverageTOEnergyThetaData","Creating %i CT bins from data 1 (%i) and data 2 (%i)",
         iCTbin,e1->GetNCosThetaBin(),e2->GetNCosThetaBin());

    TOEnergyThetaData* e3 = new TOEnergyThetaData(iEbin,iCTbin,"cs");

    for (Int_t i = 0; i < iEbin; i++)
    {
        if (i < e1->GetNEnergyBin() && i < e2->GetNEnergyBin())
        {
           if (e1->GetEnergy(i)!=e2->GetEnergy(i))
           {
              Error("TMTools::AverageTOEnergyThetaData","Not same energy in bin %i (%5.3f - %5.3f)",
              i,e1->GetEnergy(i),e2->GetEnergy(i));
           }
           else
           {
              e3->SetEnergy(i,e1->GetEnergy(i));
              e3->SetEnergyError(i, 0.5*TMath::Sqrt(e1->GetEnergyError(i)*e1->GetEnergyError(i)+
                                 e2->GetEnergyError(i)*e2->GetEnergyError(i)));
           }
        }
        else if ( i >= e1->GetNEnergyBin() )
        {
           e3->SetEnergy(i,e2->GetEnergy(i) );
           e3->SetEnergyError(i, e2->GetEnergyError(i));
        }
        else if ( i >= e2->GetNEnergyBin() )
        {
           e3->SetEnergy(i,e1->GetEnergy(i));
           e3->SetEnergyError(i, e1->GetEnergyError(i));
        }
        else
           Error("TMTools::AverageTOEnergyThetaData","WTF!");

        for (Int_t j = 0; j < iCTbin; j++)
        {
           if (j < e1->GetNCosThetaBin() && j < e2->GetNCosThetaBin())
           {
              if (e1->GetCosThetaBinCenter(j) != e2->GetCosThetaBinCenter(j))
              {
                 Error("TMTools::AverageTOEnergyThetaData","Not same CT value in bin %i (%5.3f - %5.3f)",
                 i,e1->GetCosThetaBinCenter(j),e2->GetCosThetaBinCenter(j));
              }
              else
              {
                 e3->SetValue(i,j,0.5*(e1->GetValue(i,j)+e2->GetValue(i,j)));
                 e3->SetError(i,j,0.5*TMath::Sqrt(e1->GetError(i,j)*e1->GetError(i,j)+e2->GetError(i,j)*e2->GetError(i,j)));
              }
           }
           else if ( j >= e1->GetNCosThetaBin() )
           {
              e3->SetValue(i,j,e2->GetValue(i,j));
              e3->SetError(i,j,e2->GetError(i,j));
           }
           else if ( j >= e2->GetNCosThetaBin() )
           {
              e3->SetValue(i,j,e1->GetValue(i,j));
              e3->SetError(i,j,e1->GetError(i,j));
           }
           else
              Error("TMTools::AverageTOEnergyThetaData","WTF!");
        } 
    }

    return e3;
}

//______________________________________________________________________________
TOEnergyThetaData* TMTools::AddTOEnergyThetaData(TOEnergyThetaData* e1, TOEnergyThetaData* e2)
{
    Int_t iEbin = e1->GetNEnergyBin();
    if (e2->GetNEnergyBin() != iEbin)
       Error("TMTools::AverageTOEnergyThetaData","Not same energy bins in data 1 (%i) and data 2 (%i)",
             e1->GetNEnergyBin(),e2->GetNEnergyBin());
    else
       Info("TMTools::AverageTOEnergyThetaData","Adding %i energy bins from data 1 (%i) and data 2 (%i)",iEbin,e1->GetNEnergyBin(),e2->GetNEnergyBin());

    Int_t iCTbin = e1->GetNCosThetaBin();
    if (e2->GetNCosThetaBin() != iCTbin)
       Error("TMTools::AverageTOEnergyThetaData","Not same CT bins in data 1 (%i) and data 2 (%i)",
             e1->GetNCosThetaBin(),e2->GetNCosThetaBin());
    else
       Info("TMTools::AverageTOEnergyThetaData","Adding %i CT bins from data 1 (%i) and data 2 (%i)",
            iCTbin,e1->GetNCosThetaBin(),e2->GetNCosThetaBin());

    TOEnergyThetaData* e3 = new TOEnergyThetaData(iEbin,iCTbin,"cs");

    for (Int_t i = 0; i < iEbin; i++)
    {
        if (e1->GetEnergy(i)!=e2->GetEnergy(i))
        {
           Error("TMTools::AverageTOEnergyThetaData","Not same energy in bin %i (%5.3f - %5.3f)",
           i,e1->GetEnergy(i),e2->GetEnergy(i));
        }
        else
        {
           e3->SetEnergy(i,e1->GetEnergy(i));
           e3->SetEnergyError(i, e1->GetEnergyError(i));
        }

        for (Int_t j = 0; j < iCTbin; j++)
        {
           if (e1->GetCosThetaBinCenter(j) != e2->GetCosThetaBinCenter(j))
           {
              Error("TMTools::AverageTOEnergyThetaData","Not same CT value in bin %i (%5.3f - %5.3f)",
              i,e1->GetCosThetaBinCenter(j),e2->GetCosThetaBinCenter(j));
           }
           else
           {
              e3->SetValue(i,j,e1->GetValue(i,j)+e2->GetValue(i,j));
              e3->SetError(i,j,TMath::Sqrt(e1->GetError(i,j)*e1->GetError(i,j)+e2->GetError(i,j)*e2->GetError(i,j)));
           }
        } 
    }

    return e3;
}

//______________________________________________________________________________
void TMTools::PrintModelData(TGraphErrors* g, const Char_t* fModel)
{
    FILE* fData;
    fData = fopen(Form("%s.dat",fModel),"w");

    fprintf(fData, "%s - W [MeV]; TCS[mub]\n",fModel);
    fprintf(fData, "\n");

    for (Int_t i = 0; i < g->GetN(); i++)
    {
       Double_t e = 0.;
       Double_t d = 0.;

       g->GetPoint(i, e, d);
       fprintf(fData, "%5.3f %5.3f\n", e, d);
    }

    fclose(fData);

    return;
}

//______________________________________________________________________________
void TMTools::PrintData(TGraphErrors* g, const Char_t* fModel)
{
    FILE* fData;
    fData = fopen(Form("%s.dat",fModel),"w");

    fprintf(fData, "%s - W [MeV]; TCS[mub]; eTCS[mub]\n",fModel);
    fprintf(fData, "\n");

    for (Int_t i = 0; i < g->GetN(); i++)
    {
       Double_t e = 0.;
       Double_t d = 0.;

       g->GetPoint(i, e, d);
       Double_t er = g->GetErrorY(i);
       fprintf(fData, "%5.3f %5.3f %5.3f\n", e, d, er);
    }

    fclose(fData);

    return;
}

//______________________________________________________________________________
Bool_t TMTools::IsObjectType(TObject* o, const Char_t* szType)
{
   if (!strcmp(o->ClassName(),szType))
      return kTRUE;
   else
      return kFALSE;
}

//______________________________________________________________________________
Bool_t TMTools::CheckObject(const Char_t* s, TFile* f)
{
   Bool_t isInFile = kFALSE;

   TIter next(f->GetListOfKeys());
   TKey *key;
   while ((key=(TKey*)next())) {
      if (!strcmp(key->GetName(),s))
         isInFile = kTRUE;
   }

   return isInFile;
}

////______________________________________________________________________________
//TGraphErrors** TMTools::GetLegendreRatios(TOEnergyThetaData* e, TList& list, Int_t pol, 
//                                          Int_t eStart, Int_t eEnd, Double_t left, Double_t right)
//{
//    const Int_t nPol = pol+1;
//    TGraphErrors* g[nPol];
//
//    TList* l = e->FitLegendre(eStart, eEnd, pol, &list, 0, left, right);
//    for (Int_t i = 0; i < nPol; i++)
//       g[i] = (TGraphErrors*) l->At(i);
//    TGraphErrors* r[pol];
//    for (Int_t i = 0; i < pol; i++)
//       r[i] = (TGraphErrors*) TOHUtils::DivideGraphs(g[i+1], g[0]);
//
//    return &r;
//}

//______________________________________________________________________________
