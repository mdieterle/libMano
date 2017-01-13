/*************************************************************************
 * Author: Manuel Dieterle, 2012
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMPluto                                                              //
//                                                                      //
// Class for easy evgen processing with Pluto++                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TMPluto.h"

ClassImp(TMPluto)

//______________________________________________________________________________
TMPluto::TMPluto(TString reac)
{
   // pass globals
   szReac = reac;
   szDis  = "1.";
   energyMin = 0.4;
   energyMax = 1.6;
   nEvents = 100000;
   nFiles = 1;//nfi;
   fName = NULL;//filename;
   fLoc = "./";//filelocation;
   IsMkin = kTRUE;
   fDiam = 1.0;
   fLength = 1.0;
   IsMkin = kFALSE;

   // init main particles
   strBeam = "";
   strTarget = "";
   strPart = "";
   strSpec = "";

   // init generation counter
   nGen = 0;
   iIter = 0;

   // init particles per generation counter
   nPar = new Int_t[10];
   // init particles per decay counter
   iHierarchy = new Int_t[10];

   // init particle per generation array
   strDec = new TString*[10];
   // init particle per decay array
   strHierarchy = new TString*[10];

   // init generation dependance
   iOrigGen = new Int_t*[10];
   // init particle dependence
   iOrigPar = new Int_t*[10];

   for (Int_t i = 0; i < 10; i++)
   {
      nPar[i] = 0;
      iHierarchy[i] = 0;
      strDec[i] = new TString[50];
      strHierarchy[i] = new TString[50];
      iOrigGen[i] = new Int_t[50];
      iOrigPar[i] = new Int_t[50];

      for (Int_t j = 0; j < 50; j++)
      {
         strDec[i][j] = "";
         strHierarchy[i][j] = "";
         iOrigGen[i][j] = 0;
         iOrigPar[i][j] = 0;
      }
   }
   
   // set initial state particles
   this->SetIS(szReac);
   // set final state particles
   this->SetFS(szReac);
   // set decay hierarchy
   this->InitHierarchy();
   // print informations
   this->Print();
   // generate pluto evgen
//   this->Generate();
}

//______________________________________________________________________________
TMPluto::~TMPluto()
{
   for (Int_t i = 0; i < 10; i++)
   {
      for (Int_t j = 0; j < 50; j++)
      {
         delete strDec[i][j];
         delete strHierarchy[i][j];
         delete &iOrigGen[i][j];
         delete &iOrigPar[i][j];
      }

      delete &nPar[i];
      delete &iHierarchy[i];
      delete [] strDec[i];
      delete [] strHierarchy[i];
      delete [] iOrigGen[i];
      delete [] iOrigPar[i];
   }

   delete [] nPar;
   delete [] iHierarchy;
   delete [] strDec;
   delete [] strHierarchy;
   delete [] iOrigGen;
   delete [] iOrigPar;
}

//______________________________________________________________________________
void TMPluto::SetIS(TString szReac)
{
     // get and remove beam particle
     strBeam = szReac(0,szReac.First('+')-1);
     szReac.Remove(0,szReac.First('+')+2);

     // get target particle
     strTarget = szReac(0,szReac.First('-')-1);

     // if target is quasi-free
     // set participant and spectator
     if (strTarget.Contains("("))
     {
        strPart   = strTarget(strTarget.First('(')+1,(strTarget.First(')')-1)-(strTarget.First('(')));
        strSpec   = strTarget(strTarget.Last('(')+1,(strTarget.Last(')')-1)-(strTarget.Last('(')));
        strTarget = strTarget(0,strTarget.First('('));
     }

   return;
}

//______________________________________________________________________________
void TMPluto::SetFS(TString szReac)
{
   // remove IS
   szReac.Remove(0,szReac.First('>')+2);

   Int_t nOrigPar = -1;
   Int_t nGenTmp = 0;

   while ( szReac != "" )
   {
      Int_t iNext = 0;

      while ( (szReac.First(' ') == 0) ||
              ((szReac.First('+') == 0) && ((szReac.First(' ') == 1))) ||
              (szReac.First('[') == 0) ||
              (szReac.First(']') == 0) )
      {
         if (szReac.First('[') == 0)
            nGenTmp++;
         else if (szReac.First(']') == 0)
            nGenTmp--;
         szReac.Remove(0,1);
      }

      Int_t iPlus  = szReac.First('+');
      if (iPlus < 0) iPlus = 9999;
      if ((szReac.First('+')-1) != szReac.First(' ')) iPlus = 9999;
      if ((szReac.First('+')) < szReac.First(' ')) iPlus = szReac.First('+')+2;
//Printf("Can not be (%s): %d", szReac.Data(), szReac.First('+ +'));
      Int_t iLeft  = szReac.First('[');
      if (iLeft < 0) iLeft = 9999;
      Int_t iRight = szReac.First(']');
      if (iRight < 0) iRight = 9999;

      Bool_t kPlus  = kFALSE;
      Bool_t kLeft  = kFALSE;
      Bool_t kRight = kFALSE;
      Bool_t kEnd   = kFALSE;

      if ( (iPlus < iLeft) && 
           (iPlus < iRight) )
      {
         iNext = iPlus;
         kPlus = kTRUE;
      }
      else if ( (iLeft < iPlus) && 
                (iLeft < iRight) )
      {
         iNext = iLeft;
         kLeft = kTRUE;
      }
      else if ( (iRight < iPlus) &&
                (iRight < iLeft) )
      {
         iNext = iRight;
         kRight = kTRUE;
      }
      else if ( (iPlus  == 9999) &&
                (iLeft  == 9999) &&
                (iRight == 9999) )
      {
         iNext = szReac.Length() - 1;
         kEnd = kTRUE;
      }      
      
      if (kPlus)
      {
         strDec[nGenTmp][nPar[nGenTmp]] = szReac(0,iNext-1);
         iOrigGen[nGenTmp][nPar[nGenTmp]] = nGenTmp - 1;
         iOrigPar[nGenTmp][nPar[nGenTmp]] = nPar[nGenTmp-1]-1;
         szReac.Remove(0,iNext+2);
         nPar[nGenTmp]++;
      }
      else if (kLeft)
      {
         strDec[nGenTmp][nPar[nGenTmp]] = szReac(0,iNext);
         iOrigGen[nGenTmp][nPar[nGenTmp]] = nGenTmp - 1;
         iOrigPar[nGenTmp][nPar[nGenTmp]] = nPar[nGenTmp-1]-1;
         nOrigPar = nPar[nGenTmp];
         szReac.Remove(0,iNext+1);
         nPar[nGenTmp]++;
         nGenTmp++;
      }
      else if (kRight)
      {
         strDec[nGenTmp][nPar[nGenTmp]] = szReac(0,iNext);
         iOrigGen[nGenTmp][nPar[nGenTmp]] = nGenTmp - 1;
         iOrigPar[nGenTmp][nPar[nGenTmp]] = nPar[nGenTmp-1]-1;
         nOrigPar = nPar[nGenTmp-1];
         szReac.Remove(0,iNext+1);
         nPar[nGenTmp]++;
         nGenTmp--;
      }
      else if (kEnd)
      {
         strDec[nGenTmp][nPar[nGenTmp]] = szReac;
         iOrigGen[nGenTmp][nPar[nGenTmp]] = nGenTmp - 1;
         iOrigPar[nGenTmp][nPar[nGenTmp]] = nPar[nGenTmp-1]-1;
         szReac = "";
         nPar[nGenTmp]++;
      }
      if (nGen < nGenTmp)
         nGen = nGenTmp;
   }
}

//______________________________________________________________________________

void TMPluto::InitHierarchy()
{
   if (strPart != "")
   {
      iIter = 1;

      iHierarchy[0] = 5;
      strHierarchy[0][0] = strBeam;
      strHierarchy[0][1] = strTarget;
      strHierarchy[0][2] = strBeam;
      strHierarchy[0][3] = strPart;
      strHierarchy[0][4] = strSpec;

      iHierarchy[1] = 2;
      strHierarchy[1][0] = strBeam;
      strHierarchy[1][1] = strPart;
   }
   else
   {
      iIter = 0;

      iHierarchy[0] = 2;
      strHierarchy[0][0] = strBeam;
      strHierarchy[0][1] = strTarget;
   }

   for (Int_t j = 0; j < nPar[0]; j++)
   {
       iHierarchy[iIter] ++;
       strHierarchy[iIter][2+j] = strDec[0][j];
   } 
 
   for (Int_t i = 0; i < nGen+1; i++)
   {
      for (Int_t j = 0; j < nPar[i]; j++)
      {
         Bool_t IsDec = kFALSE;
         Int_t nDec = 0;
         for (Int_t k = 0; k < nPar[i+1]; k++)
         {
            if ( (iOrigGen[i+1][k] == i) && (iOrigPar[i+1][k] == j) )
            {
               IsDec = kTRUE;
               nDec = k;
            }
         }
         if (IsDec)
         {
            iIter++;
            iHierarchy[iIter] = 1;
            strHierarchy[iIter][0] = strDec[i][j];

            Int_t n = 1;
            for (Int_t k = 0; k < nPar[i+1]; k++)
            {
               if ( (iOrigGen[i+1][k] == i) && (iOrigPar[i+1][k] == j) )
               {
                  iHierarchy[iIter]++;
                  strHierarchy[iIter][n] = strDec[i+1][k];
                  n++;
               }
            }
         }
      }
   }

   if (strPart == "")
      iIter++;

   return;
}

//______________________________________________________________________________

void TMPluto::Print()
{
   Info("TMPluto::Print","Printing all Particle Informations");
   printf("\n");

   if (strPart != "")
      printf("1. %s + %s --> (%s + %s) + %s\n", strHierarchy[0][0].Data(),
                                               strHierarchy[0][1].Data(),
                                               strHierarchy[0][2].Data(),
                                               strHierarchy[0][3].Data(),
                                               strHierarchy[0][4].Data());

   Int_t ii = 0;
   
   for (Int_t i = 0; i < iIter; i++)
   {
      if (strPart != "")
         ii = i+1;
      else
         ii = i;
 
      if (i == 0)
      {
          printf("%i. %s + %s --> ", ii+1, strHierarchy[ii][0].Data(),
                                          strHierarchy[ii][1].Data());

          for (Int_t j = 2; j < iHierarchy[ii]-1; j++)
             printf("%s + ", strHierarchy[ii][j].Data());

          printf("%s\n", strHierarchy[ii][iHierarchy[ii]-1].Data());
      }
      else
      {
         printf("%i. %s --> ", ii+1, strHierarchy[ii][0].Data());

         for (Int_t j = 1; j < iHierarchy[ii]-1; j++)
            printf("%s + ", strHierarchy[ii][j].Data());

         printf("%s\n", strHierarchy[ii][iHierarchy[ii]-1].Data());
      }
   }
   printf("\n");
}

//______________________________________________________________________________
void TMPluto::Generate()
{
   Char_t szName[256];

   if (strPart != "")
   {
      makeDistributionManager()->Exec("nucleus_fermi:gamma");
      Info("TMPluto::Generate","Enabling Quasi-Free Production");
   }

   //call member of the modified PBeamSmearing class
   PBeamSmearing *smear = new PBeamSmearing("beam_smear", "Beam Momentum Smearing");

   sprintf(szName,"%s + %s",strBeam.Data(),strTarget.Data());

   //pass reaction
   smear->SetReaction(szName);
   //define function + momentum range according to which the beam momentum is smeared
   smear->SetMomentumFunction(new TF1(szDis,szDis,energyMin,energyMax));
   //add smearing to pluto world
   makeDistributionManager()->Add(smear);

   for (Int_t k = 0; k < nFiles; k++)
   {
      PUtils::SetSeed(0);

      PParticle *pBeam   = new PParticle(strBeam.Data(),energyMin);
      PParticle *pTarget = new PParticle(strTarget.Data());
      PParticle *s       = new PParticle(*pBeam+*pTarget);

      PParticle ***cc;
      cc = new PParticle**[iIter+1];

      PChannel** c;
      c = new PChannel*[iIter+1];

      for (Int_t i = 0; i < iIter+1; i++)
      {
         c[i] = 0;
         cc[i] = new PParticle*[10];

         for (Int_t j = 0 ; j < 10; j++)
         {
            cc[i][j] = 0;
         }
      }

      PParticle *pBeam2 = 0;
      PParticle *pParticipant = 0;
      PParticle *s2 = 0;
      PParticle *pSpectator = 0;

      Int_t iIndex = 0;
      Int_t iStart = 0;

      if (strPart != "")
      {
         pBeam2       = new PParticle(strBeam.Data());
         pParticipant = new PParticle(strPart.Data());
         s2           = new PParticle(*pBeam2+*pParticipant);

         pSpectator = new PParticle(strSpec.Data());

         cc[0][0] = s;
         cc[0][1] = pSpectator;
         cc[0][2] = s2;
         cc[1][0] = s2;

         c[0] = new PChannel(cc[0],2);

         iIndex = 1;
         iStart = 1;
      }
      else
      {
         cc[0][0] = s;
      }

      PParticle* pTmp = 0;

      for (Int_t i = 0; i < iIter; i++)
      {
         if (strPart != "")
            iStart = i+1;
         else
            iStart = i;

         if (i == 0)
         {
            for (Int_t j = 2; j < iHierarchy[iStart]; j++)
            {
               pTmp = new PParticle(strHierarchy[iStart][j]);
               cc[iIndex][j-1] = pTmp;
            }
         }
         else
         {
            for (Int_t j = 0; j < iHierarchy[iStart]; j++)
            {
               pTmp = new PParticle(strHierarchy[iStart][j]);
               cc[iIndex][j] = pTmp;
            }
         }

         if (i == 0)
            c[iIndex] = new PChannel(cc[iIndex],iHierarchy[iStart]-2);
         else
            c[iIndex] = new PChannel(cc[iIndex],iHierarchy[iStart]-1);

         iIndex++;
      }

      if (fName == NULL) fName="tmp";
      PReaction *Reac = new PReaction(c,Form("%s_%i",fName,k),iStart+1,1,0,0,0);    // Define reaction
///   
      Reac->Print();        // Write to .root file 
      Reac->loop(nEvents);

      if (IsMkin)
      {
         // start pluto2mkin for 1.76cm target
         gSystem->Exec(Form("pluto2mkin --beam diam=%5.3f --target length=%5.3f --input %s_%i.root",fDiam,fLength,fName,k));
   
         // delete pluto file
         gSystem->Unlink(Form("%s_%i.root",fName,k));
   
         // rename mkin file
         gSystem->Rename(Form("%s_%i_mkin.root",fName,k), Form("%s_%i.root",fName,k));
      }
 
      // move to desired location
      gSystem->Exec(Form("mv %s_%i.root %s",fName,k,fLoc));
   
      printf("\n File %s_%i.root moved to %s\n",fName,k,fLoc); 

      delete Reac;
//      for (Int_t i = 0; i < iIter+1; i++)
//      {
//         for (Int_t j = 0; j < 10; j++)
//         {
//            delete cc[i][j];
//         }
//
//         delete c[i];
//         delete [] cc[i];
//      }
//
//      delete [] c;
//      delete [] cc;
      delete pBeam;
      delete pTarget;
      delete s;
      delete pTmp;

      if (strPart != "")
      {
         delete pBeam2;
         delete pParticipant;
         delete pSpectator;
         delete s2;
      }
   }
}

//______________________________________________________________________________


