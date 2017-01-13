/*************************************************************************
 * Author: Manuel Dieterle, 2012
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMObjectCollection                                                   //
//                                                                      //
// Class for collecting any kind of TObjects                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TMObjectCollection.h"

ClassImp(TMObjectCollection)

//______________________________________________________________________________
TMObjectCollection::TMObjectCollection(const Char_t* szName, const Char_t* szTitle)
    : TMCollection(szName, szTitle)
{
    // Constructor. Allocate space for 'nMaxTaggCh' tagger channels.
 
    // init members
    strDat   = "";
    strEmpty = "";
    strSig   = "";
    strSig1  = "";
    strSig2  = "";
    strBG    = "";
    strTot   = "";
    strBin   = "";
    strArr   = "";
    strFit   = "";

    nArr = 0;
    nCT = 0;
    iCT = 0;
    nMM = 0;
    iMM = 0;

    yMin = 0.;
    yMax = 0.;
}

//______________________________________________________________________________
TMObjectCollection::TMObjectCollection()
    : TMCollection()
{
}
//______________________________________________________________________________
TMObjectCollection::~TMObjectCollection()
{
    // Destructor.
//    if (hArrHisto) delete [] hArrHisto;
}

//______________________________________________________________________________
void TMObjectCollection::AddEnergyBin(TMEBin* e)
{
    strBin = e->GetName();
    fCollection->Add(e);    

    return;
}

//______________________________________________________________________________
void TMObjectCollection::Add3DBin(TMEBin* e, Int_t nCTbin, Int_t iCTbin, Int_t nMMbin, Int_t iMMbin)
{
    if (strBin != "")
       Error("TMObjectCollection::Add3DBin","Already added TMEBin %s\n",strBin.Data());
    else
    {
       strBin = e->GetName();
       fCollection->Add(e);
       nCT = nCTbin;
       iCT = iCTbin;
       nMM = nMMbin;
       iMM = iMMbin;
    }

    return;
}

//______________________________________________________________________________
void TMObjectCollection::AddDataObject(TObject* o)
{
    // Set Signal Histogram
    strDat = o->GetName();
    fCollection->Add(o);

    return;
}

//______________________________________________________________________________
void TMObjectCollection::AddEmptyTargetObject(TObject* o)
{
    // Set Signal Histogram
    strEmpty = o->GetName();
    fCollection->Add(o);

    return;
}

//______________________________________________________________________________
void TMObjectCollection::AddSignalObject(TObject* o)
{
    // Set Signal Histogram
    strSig = o->GetName();
    fCollection->Add(o);

    return;
}

//______________________________________________________________________________
void TMObjectCollection::AddSignal1Object(TObject* o)
{
    // Set Signal Histogram
    strSig1 = o->GetName();
    fCollection->Add(o);

    return;
}

//______________________________________________________________________________
void TMObjectCollection::AddSignal2Object(TObject* o)
{
    // Set Signal Histogram
    strSig2 = o->GetName();
    fCollection->Add(o);

    return;
}

//______________________________________________________________________________
void TMObjectCollection::AddBackgroundObject(TObject* o)
{
    // Set Signal Histogram
    strBG = o->GetName();
    fCollection->Add(o);

    return;
}

//______________________________________________________________________________
void TMObjectCollection::AddTotalObject(TObject* o)
{
    // Set Signal Histogram
    strTot = o->GetName();
    fCollection->Add(o);

    return;
}

//______________________________________________________________________________
void TMObjectCollection::AddFitObject(TObject* o)
{
    // Set Signal Histogram
    strFit = o->GetName();
    fCollection->Add(o);

    return;
}

//______________________________________________________________________________
void TMObjectCollection::AddArrayObject(Int_t n, TObject** o)
{
    // Set Signal Histogram
    nArr = n;
    Char_t strArr_tmp[256];
    if (TMTools::IsObjectType(o[0],"TFile*"))
       sscanf(o[0]->GetName(),"%[^_]_%*i.root",strArr_tmp);
    else
       sscanf(o[0]->GetName(),"%[^_]_%*i",strArr_tmp);
    strArr = strArr_tmp;
    for (Int_t i = 0; i < n; i++)
       fCollection->Add(o[i]);
    return;
}

//______________________________________________________________________________
TObject** TMObjectCollection::GetArrayObject()
{
   TObject** o;
   o = new TObject*[nArr];
   Char_t s[256];

   for (Int_t i = 0; i < nArr; i++)
   {
      sprintf(s,"%s_%i",strArr.Data(),i);

      if(!(o[i]=(TObject*)GetObject(s)))
      {
         sprintf(s,"%s_%i.root",strArr.Data(),i);
         o[i]=(TObject*)GetObject(s);
      }
   }
   
   return o;
}
//______________________________________________________________________________
TObject* TMObjectCollection::GetArrayObject(Int_t i)
{
   TObject* o;
   Char_t s[256];
   Int_t j;
   Int_t k;

   sscanf(this->GetName(),"%*[^_]_%i_%i",&j,&k);
   sprintf(s,"%s_%i_%i_%i",strArr.Data(),i,j,k);
  
   if(!(o=((TObject*)GetObject(s))))
   {
      sprintf(s,"%s%i_%i_%i.root",strArr.Data(),i,j,k);
      o=((TObject*)GetObject(s));
      printf("  Exporting %s\n",s);
   }

   return o;
}

//______________________________________________________________________________
void TMObjectCollection::Print()
{
    // Print out the content of this class.
    this->ListObjects();
//    printf("Found Histos in this collection:\n");
//    if (strDat!="") printf("    Data:    %s\n", strDat.Data());
//    if (strSig1!="") printf("    Signal1: %s\n", strSig1.Data());
//    if (strSig2!="") printf("    Signal2: %s\n", strSig2.Data());
//    if (strBG!="") printf("    BG:      %s\n", strBG.Data());
//    if (strTot!="") printf("    Total:   %s\n", strTot.Data());
//    if (strBGf!="")
//    {
//       printf("Found BG Function in this collection:\n");
//       printf("    BG: %s\n", strBGf.Data());
//    }
    if (strBin!="")
    {
       printf("Energy                      : %f\n", ((TMEBin*)GetObject(strBin.Data()))->GetEnergy());
       printf("Energy error                : %f\n", ((TMEBin*)GetObject(strBin.Data()))->GetEnergyError());
       printf("Total # of tagg. ch.        : %d\n", ((TMEBin*)GetObject(strBin.Data()))->GetNTotalTaggerChannels());
       printf("Tagger Channel bin          : %d\n", ((TMEBin*)GetObject(strBin.Data()))->GetNTaggerChannel());
       printf("# of tagg. ch. in this bin  : %d\n", ((TMEBin*)GetObject(strBin.Data()))->GetNTaggerChannels());
       printf("first tagg. ch. in this bin : %d\n", ((TMEBin*)GetObject(strBin.Data()))->GetTaggerChannelLo());
       printf("last tagg. ch. in this bin  : %d\n", ((TMEBin*)GetObject(strBin.Data()))->GetTaggerChannelUp());
    }
    if (nCT)
       printf("Total # of CT bins          : %d\n", nCT);
    if (iCT)
       printf("Cos(Theta) bin              : %d\n", iCT);
    if (nMM)
       printf("Total # of MM bins          : %d\n", nMM);
    if (iMM)
       printf("Missing Mass bin            : %d\n", iMM);
    return;
}

//______________________________________________________________________________
void TMObjectCollection::GetHistoRange()
{
   if (strcmp(strTot,"") && strcmp(strDat,""))
   {
      if (!strcmp(GetObject(strTot.Data())->ClassName(),"TF1"))
      {
         yMin = ((TH1*)GetObject(strDat.Data()))->GetMinimum();
         yMax = ((TH1*)GetObject(strDat.Data()))->GetMaximum();
         if ((((TF1*)GetObject(strTot.Data()))->GetMinimum())<yMin)
            yMin = ((TF1*)GetObject(strTot.Data()))->GetMinimum();
         if ((((TF1*)GetObject(strTot.Data()))->GetMaximum())>yMax)
            yMax = ((TF1*)GetObject(strTot.Data()))->GetMaximum();
      }
      else
      {
         yMin = ((TH1*)GetObject(strDat.Data()))->GetMinimum();
         yMax = ((TH1*)GetObject(strDat.Data()))->GetMaximum();
         if ((((TH1*)GetObject(strTot.Data()))->GetMinimum())<yMin)
            yMin = ((TH1*)GetObject(strTot.Data()))->GetMinimum();
         if ((((TH1*)GetObject(strTot.Data()))->GetMaximum())>yMax)
            yMax = ((TH1*)GetObject(strTot.Data()))->GetMaximum();
      }
   }
   else if (strcmp(strTot,"") && !strcmp(strDat,""))
   {
      Error("TMObjectCollection::GetHistoRange","Data Histogram not found");
   }                                       
   else if (!strcmp(strTot,"") && strcmp(strDat,""))
   {
      Error("TMObjectCollection::GetHistoRange","Total Histogram not found");
   }
   else
      Error("TMObjectCollection::GetHistoRange","Data and Total Histogram not found");

   yMin *= 1.2;
   yMax *= 1.2;

   return;
}

//______________________________________________________________________________
void TMObjectCollection::Draw(const Char_t* szType, Double_t xMin, Double_t xMax)
{
   if (strcmp(szType,"all")) 
   {
      Info("TMObjectCollection::Draw","Not all necessary histograms available");
   }
   else
   {
      this->GetHistoRange();
      if (!strcmp(GetObject(strBG.Data())->ClassName(),"TF1"))
      {
         ((TF1*)GetObject(strBG.Data()))->Draw();
         ((TF1*)GetObject(strBG.Data()))->GetXaxis()->SetRangeUser(xMin,xMax);
         ((TF1*)GetObject(strBG.Data()))->GetYaxis()->SetRangeUser(yMin,yMax);
         ((TH1*)GetObject(strSig.Data()))->Draw("Hsame");
         ((TH1*)GetObject(strTot.Data()))->Draw("Hsame");
         ((TH1*)GetObject(strDat.Data()))->SetLineColor(1);
         ((TH1*)GetObject(strDat.Data()))->Draw("E1same");
      }
      else
      {
         ((TH1*)GetObject(strBG.Data()))->Draw("H");
         ((TH1*)GetObject(strBG.Data()))->GetXaxis()->SetRangeUser(xMin,xMax);
         ((TH1*)GetObject(strBG.Data()))->GetYaxis()->SetRangeUser(yMin,yMax);
         ((TH1*)GetObject(strSig.Data()))->Draw("Hsame");
         ((TH1*)GetObject(strTot.Data()))->Draw("Hsame");
         ((TH1*)GetObject(strDat.Data()))->SetLineColor(1);
         ((TH1*)GetObject(strDat.Data()))->Draw("E1same");
      }
   }
   return;
}

//______________________________________________________________________________
void TMObjectCollection::DrawPaperCollection(Double_t xMin, Double_t xMax)
{
   this->GetHistoRange();

   // Draw Data Histo
   ((TH1*)this->GetDataObject())->SetLineColor(1);
   ((TH1*)this->GetDataObject())->SetTitle("");
   ((TH1*)this->GetDataObject())->SetStats(0);
   ((TH1*)this->GetDataObject())->Draw("E1");
   ((TH1*)this->GetDataObject())->GetXaxis()->SetRangeUser(xMin,xMax);
   ((TH1*)this->GetDataObject())->GetYaxis()->SetRangeUser(yMin,yMax);

   // Draw Total Background Object
   if (strcmp(strBG,""))
   {
      if (!strcmp(GetObject(strTot.Data())->ClassName(),"TF1"))
      {
         ((TF1*)this->GetBackgroundObject())->SetLineColor(3);
         ((TF1*)this->GetBackgroundObject())->Draw("Hsame");
      }
      else
      {
         ((TH1*)this->GetBackgroundObject())->SetLineColor(3);
         ((TH1*)this->GetBackgroundObject())->Draw("Hsame");
      }
   }
   
   // Draw Background Contributions
   if (strcmp(strArr,""))
   {
      for (Int_t i = 0; i < nArr; i++)
      {
         ((TH1*)this->GetArrayObject(i))->SetLineColor(7);
         ((TH1*)this->GetArrayObject(i))->Draw("Hsame");
      }      
   }

   // Draw Signal Object
   if (strcmp(strSig,""))
   {
      ((TH1*)this->GetSignalObject())->SetLineColor(4);
      ((TH1*)this->GetSignalObject())->Draw("Hsame");
   }
//   else if (strcmp(strSig1,"") && strcmp(strSig2,""))
//   {
//      ((TH1*)this->GetSignal1Object())->SetLineColor(4);
//      ((TH1*)this->GetSignal1Object())->Draw("Hsame");
//      ((TH1*)this->GetSignal2Object())->SetLineColor(6);
//      ((TH1*)this->GetSignal2Object())->Draw("Hsame");
//   }
//   else
//      Error("TMObjectCollection::DrawFitCollection","No Signal Object found in Collection");

   // Draw Total Fit Object
   if (strcmp(strTot,""))
   {
      ((TH1*)this->GetTotalObject())->SetLineColor(2);
      ((TH1*)this->GetTotalObject())->Draw("Hsame");
   }

   // Draw Fit Function
   if (strcmp(strFit,""))
   {
      ((TF1*)this->GetFitObject())->SetLineColor(6);
      ((TF1*)this->GetFitObject())->SetLineStyle(2);
      ((TF1*)this->GetFitObject())->Draw("same");
   }

   Char_t szE[256];
   sprintf(szE,"#font[41]{%s}",(Char_t*)this->GetTitle());

   TLatex latex;
   latex.SetTextSize(0.1);
   latex.DrawLatex(xMin+0.1*(xMax-xMin),0.9*yMax,szE);
}



//______________________________________________________________________________
void TMObjectCollection::DrawFitCollection(Double_t xMin, Double_t xMax)
{
   this->GetHistoRange();

   ((TH1*)this->GetDataObject())->SetLineColor(1);
   ((TH1*)this->GetDataObject())->SetTitle("");
   ((TH1*)this->GetDataObject())->SetStats(0);
   ((TH1*)this->GetDataObject())->Draw("E1");
   ((TH1*)this->GetDataObject())->GetXaxis()->SetRangeUser(xMin,xMax);
   ((TH1*)this->GetDataObject())->GetYaxis()->SetRangeUser(yMin,yMax);

   if (strcmp(strBG,""))
   {
      if (!strcmp(GetObject(strTot.Data())->ClassName(),"TF1"))
      {
         ((TF1*)this->GetBackgroundObject())->SetLineColor(7);
         ((TF1*)this->GetBackgroundObject())->Draw("Hsame");
      }
      else
      {
         ((TH1*)this->GetBackgroundObject())->SetLineColor(7);
         ((TH1*)this->GetBackgroundObject())->Draw("Hsame");
      }
   }
   
   if (strcmp(strArr,""))
   {
//printf("Number of background contributions %i\n",nArr);
      for (Int_t i = 0; i < nArr; i++)
      {
         if (i != nArr-1)
         {
            ((TH1*)this->GetArrayObject(i))->SetLineColor(3);
            ((TH1*)this->GetArrayObject(i))->Draw("Hsame");
         }
         else
         {
            ((TH1*)this->GetArrayObject(i))->SetLineColor(3);
//            ((TH1*)this->GetArrayObject(i))->SetLineColor(7);
            ((TH1*)this->GetArrayObject(i))->Draw("Hsame");
         }
      }      
   }

   if (strcmp(strSig,""))
   {
      ((TH1*)this->GetSignalObject())->SetLineColor(4);
      ((TH1*)this->GetSignalObject())->Draw("Hsame");
   }
   else if (strcmp(strSig1,"") && strcmp(strSig2,""))
   {
      ((TH1*)this->GetSignal1Object())->SetLineColor(4);
      ((TH1*)this->GetSignal1Object())->Draw("Hsame");
      ((TH1*)this->GetSignal2Object())->SetLineColor(6);
      ((TH1*)this->GetSignal2Object())->Draw("Hsame");
   }
//   else
//      Error("TMObjectCollection::DrawFitCollection","No Signal Object found in Collection");
   ((TH1*)this->GetTotalObject())->SetLineColor(2);
   ((TH1*)this->GetTotalObject())->Draw("Hsame");
   if (strcmp(strFit,""))
   {
      ((TF1*)this->GetFitObject())->SetLineColor(4);
      ((TF1*)this->GetFitObject())->SetLineStyle(2);
      ((TF1*)this->GetFitObject())->Draw("same");
   }

   Char_t szE[256];
   sprintf(szE,"#font[41]{%s}",(Char_t*)this->GetTitle());

   TLatex latex;
   latex.SetTextSize(0.1);
   latex.DrawLatex(xMin+0.1*(xMax-xMin),0.9*yMax,szE);
}

//______________________________________________________________________________
void TMObjectCollection::DrawArrayCollection(Int_t iRebin, Double_t xMin, Double_t xMax)
{
   Double_t yMin = 0.;
   Double_t yMax = 0.;

   for (Int_t i = 0; i < this->GetArraySize(); i++)
   {
      ((TH1*)this->GetArrayObject(i))->Rebin(iRebin);

      if (((TH1*)this->GetArrayObject(i))->GetMinimum() < yMin) 
         yMin = ((TH1*)this->GetArrayObject(i))->GetMinimum();
      if (((TH1*)this->GetArrayObject(i))->GetMaximum() > yMax) 
         yMax = ((TH1*)this->GetArrayObject(i))->GetMaximum();
   }

   for (Int_t i = 0; i < this->GetArraySize(); i++)
   {
      if (i == 0)
      {
         ((TH1*)this->GetArrayObject(i))->GetXaxis()->SetRangeUser(xMin,xMax);
         ((TH1*)this->GetArrayObject(i))->GetYaxis()->SetRangeUser(yMin*1.2,yMax*1.2);
         ((TH1*)this->GetArrayObject(i))->Draw("E1");
      }
      ((TH1*)this->GetArrayObject(i))->Draw("Hsame");
   }
}
//______________________________________________________________________________
