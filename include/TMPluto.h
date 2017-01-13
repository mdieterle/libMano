
/*************************************************************************
 * Author: Manuel Dieterle, 2012
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMPluto                                                         //
//                                                                      //
// Class for a collection of TObjects                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef myTMPluto
#define myTMPluto

#include "TNamed.h"
#include "PParticle.h"
#include "PDistributionManager.h"
#include "PBeamSmearing.h"
#include "PChannel.h"
#include "PReaction.h"
#include "TSystem.h"

class TMPluto : public TNamed
{

protected:

    TString szReac;
    const Char_t* szDis;
    Double_t energyMin;
    Double_t energyMax;
    Int_t nEvents;
    Int_t nFiles;
    const Char_t* fName;
    const Char_t* fLoc;
    Double_t fDiam;
    Double_t fLength;
    Bool_t IsMkin;

    TString strBeam;
    TString strTarget;
    TString strPart;
    TString strSpec;

    Int_t nGen;
    Int_t iIter;

    Int_t* nPar;
    Int_t* iHierarchy;

    TString** strDec;
    TString** strHierarchy;

    Int_t** iOrigGen;
    Int_t** iOrigPar;

public:

//    TMPluto(TString reac, const Char_t* dis, Double_t eMin, Double_t eMax, Int_t nEv);
//            Int_t nfi, const Char_t* filename, const Char_t* filelocation, 
//            Double_t diam = 1.0, Double_t length = 1.0, Bool_t convert = kTRUE);
    TMPluto(TString reac);
    virtual ~TMPluto();

    void SetBeamDistribution(const Char_t* dis){ szDis = dis; }
    void SetBeamEnergy(Double_t eMin, Double_t eMax){ energyMin = eMin; energyMax = eMax; }
    void SetNumberOfEvents(Int_t nev){ nEvents = nev; }
    void SetNumberOfFiles(Int_t nfi){ nFiles = nfi; }
    void SetFileName(const Char_t* filename){ fName = filename; }
    void SetOutputPath(const Char_t* filelocation){ fLoc = filelocation; }
    void SetBeamDiameter(Double_t diam){ fDiam = diam; }
    void SetTargetLength(Double_t length){ fLength = length; }
    void SetPluto2Mkin(Bool_t IsConv){ IsMkin = IsConv; }
    void SetIS(TString szReac);
    void SetFS(TString szReac);
    void SetMembers();
    void InitHierarchy();
    void Print();
    void Generate();

    ClassDef(TMPluto, 0)  // pluto class
};

#endif

