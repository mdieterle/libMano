/*************************************************************************
 * Author: Dominik Werthmueller, 2008-2009
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TOTools                                                               //
//                                                                      //
// Class representing an energy bin.                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TOTools.h"
#include "TOHistoCollection.h"

ClassImp(TOTools)

//______________________________________________________________________________
TH1* TOTools::Get1DHisto(Char_t* szFile, TString strOldName, TString strNewName)
{
  // get a histo "strOldName" from a file "strFile"
  // and get a clone of it "strNewName"

  TH1* h1;

  TFile *ffile = TFile::Open(szFile);
  if( ffile->IsZombie() )
  {
     cerr << " ERROR: can't open file " << ffile->GetName() << endl;
     gSystem->Exit(0);
  }

  gROOT->cd();

  h1 = (TH1*) ffile->Get( strOldName )->Clone( strNewName );

  ffile->Close();

  return h1;
}
//______________________________________________________________________________
void TOTools::DrawHistoCollection(TFile* fIn, Int_t iNPadX, Int_t iNPadY, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax, Int_t iNDC)
{
   Int_t iNPad = iNPadX*iNPadY;
   TCanvas* c = 0;
   Char_t szName[256];
   Char_t szzName[256];

   TIter next(fIn->GetListOfKeys());
   TKey *key;
   while ((key=(TKey*)next()))
   {
      if (strcmp(key->GetClassName(),"TOHistoCollection") != 0)
         printf("No object of class TOHistoCollection found\n");
      else
      {
//         printf("Found object of class TOHistoCollection of name %s\n", key->GetName());
         ((TOHistoCollection*)f->Get(key->GetName())).Print();
      }
   }




   for (Int_t i = 0; i < iCT; i++)
   {
      p[i] = new TPad*[iRTC];

      for (Int_t j = 0; j < iRTC; j++)
      {

         if (j%iNPad == 0)
         {
            if (j+iNPad>iRTC)
            {
               sprintf(szName, "Neutral_MM_CT Bin %i - Tagger Channels %03i to %03i", i, iTClo[j], iTCup[iRTC-1]);
               sprintf(szzName,"Neutral_MM_CT_Bin_%i_Tagger_Channels_%03i_to_%03i.eps", i, iTClo[j], iTCup[iRTC-1]);
            }
            else
            {
               sprintf(szName, "Neutral_MM_CT Bin %i - Tagger Channels %03i to %03i", i, iTClo[j], iTCup[j+iNPad-1]);
               sprintf(szzName,"Neutral_MM_CT_Bin_%i_Tagger_Channels_%03i_to_%03i.eps", i, iTClo[j], iTCup[j+iNPad-1]);
            }
            c = new TCanvas(szName,szName,1000,1000);
         }

         c->cd();

         sprintf(szName,"MM_CT Bin %i - Tagger Channels %03i to %03i", i, iTClo[j], iTCup[j]);

         p[i][j] = new TPad(szName, szName, 
                         xmin+(j%iNPadX)*(xmax-xmin)/iNPadX,
                         ymax-((j%iNPad)/iNPadY+1)*(ymax-ymin)/iNPadY,
                         xmin+(j%iNPadX+1)*(xmax-xmin)/iNPadX,
                         ymax-((j%iNPad)/iNPadY)*(ymax-ymin)/iNPadY);

         p[i][j]->SetLeftMargin(0);
         p[i][j]->SetRightMargin(0);
         p[i][j]->SetTopMargin(0);
         p[i][j]->SetBottomMargin(0);

         p[i][j]->Draw();
         
         p[i][j]->cd();

         hIn[i][j][0]->GetYaxis()->SetRangeUser(0.,hIn[i][j][0]->GetMaximum()*1.3);
         hIn[i][j][0]->Draw("AE1");
         for (Int_t k = 1; k < ihNum; k++)
            hIn[i][j][k]->Draw("same");

         hpAll[i][j]->Draw("same");

         TLatex aText;
         sprintf(szName, "#font[42]{E_{#gamma}=%5.3f MeV}", fMeanEnergy[j]);
         aText.SetTextSize(0.1);
         aText.DrawLatex(-850., 0.85*hIn[i][j][0]->GetMaximum(), szName);

         if (j%iNPad == iNPad-1 || j == iRTC-1)
         {
            c->cd();

            TGaxis* aX = 0;

            for (Int_t k = 0; k < iNPadX; k++)
            {
               aX = new TGaxis(xmin+k*(xmax-xmin)/iNPadX,
                                       ymax-((j%iNPad)/iNPadY+1)*(ymax-ymin)/iNPadY,
                                       xmin+(k+1)*(xmax-xmin)/iNPadX,
                                       ymax-((j%iNPad)/iNPadY+1)*(ymax-ymin)/iNPadY,
                                       ((TH1F*)hpDat[i][j])->GetXaxis()->GetBinLowEdge(1)*0.95,
                                       ((TH1F*)hpDat[i][j])->GetXaxis()->GetBinUpEdge(((TH1F*)hpDat[i][j])->GetNbinsX())*0.95,
                                       iNDC,"-");
               aX->SetLabelSize(.03);
               aX->SetLabelOffset(-.03);
               aX->SetLabelFont(42);
               aX->Draw();
            }

            TLine aLine;
            aLine.DrawLineNDC(xmax,ymax-((j%iNPad)/iNPadY+1)*(ymax-ymin)/iNPadY,xmax,ymax);

            TLatex title;
            title.SetTextAlign(22);
            title.SetTextSize(0.05);
            title.DrawLatex(xmin+0.5*(xmax-xmin),0.5*(ymin-aX->GetLabelSize())+(iNPadY-((j%iNPad)/iNPadY+1))*(ymax-ymin)/iNPadY,"#font[41]{M_{X}-M_{N} [MeV]}");
            title.SetTextSize(0.05);
            title.SetTextAngle(90);
            title.DrawLatex(0.5*xmin,ymax-0.5*((j%iNPad)/iNPadY+1)/iNPadY*(ymax-ymin),"#font[41]{Counts [a.u.]}");

            if (PSsave) c->Print(szzName);
//            fOut->cd();
            c->Write();
         }
      }
//      hExFuncD[i]->Write();
//      hExFuncP[i]->Write();
   }
//   fOut->Close();

   return;
}


}
