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


#ifndef myTMObjectCollection
#define myTMObjectCollection

#include "TMCollection.h"
#include "TNamed.h"
#include "TH1.h"
#include "TMEBin.h"
#include "TCanvas.h"
#include "TMTools.h"

class TMObjectCollection : public TMCollection
{

protected:
    Int_t nHisto;                           // number of histos
    TString strDat;                         // data object
    TString strEmpty;                       // empty target object
    TString strSig;                         // sig object
    TString strSig1;                        // sig1 object
    TString strSig2;                        // sig2 object
    TString strBG;                          // bg object
    TString strTot;                         // total object
    TString strArr;                         // array object
    TString strBin;                         // Energy Bin
    TString strFit;                         // fit func
    Int_t nArr;                             // size of array
    Int_t nCT;                              // # CTbins
    Int_t iCT;                              // CTbin
    Int_t nMM;                              // # MMbins
    Int_t iMM;                              // MMbin
    Double_t yMin;                          // ymin
    Double_t yMax;                          // ymax

public:
//    TMObjectCollection() : TMCollection() { }
//    TMObjectCollection() : TMCollection(), nHisto(0), strDat(0), strSig1(0), 
//                          strSig2(0), strBG(0), strBGf(0), strTot(0), strBin(0), 
//                          nCT(0), iCT(0), nMM(0), iMM(0), yMin(0.), yMax(0.) { }
    TMObjectCollection();
    TMObjectCollection(const Char_t* szName, const Char_t* szTitle);
    virtual ~TMObjectCollection();
    void Print();

    TObject* GetDataObject() { return (TObject*)GetObject(strDat.Data()); }
    TObject* GetEmptyTargetObject() { return (TObject*)GetObject(strEmpty.Data()); }
    TObject* GetSignalObject() { return (TObject*)GetObject(strSig.Data()); }
    TObject* GetSignal1Object() { return (TObject*)GetObject(strSig1.Data()); }
    TObject* GetSignal2Object() { return (TObject*)GetObject(strSig2.Data()); }
    TObject* GetBackgroundObject() { return (TObject*)GetObject(strBG.Data()); }
    TObject* GetTotalObject() { return (TObject*)GetObject(strTot.Data()); }
    TObject* GetFitObject() { return (TObject*)GetObject(strFit.Data()); }
//    TObject** GetArrayObject() { return (TObject**)GetObject(strArr.Data()); }
    TObject** GetArrayObject();
    TObject* GetArrayObject(Int_t i);
    Int_t GetArraySize() const { return nArr; }

    TMEBin* GetEBin() { return (TMEBin*)GetObject(strBin.Data()); }
    Int_t GetNTCBins() { return (((TMEBin*)GetObject(strBin.Data()))->GetNTotalTaggerChannels()); }
    Int_t GetTCBin() { return (((TMEBin*)GetObject(strBin.Data()))->GetNTaggerChannel()); }
    Int_t GetNCTBins() const { return nCT; }
    Int_t GetNMMBins() const { return nMM; }
    Int_t GetCTBin() const { return iCT; }
    Int_t GetMMBin() const { return iMM; }

    void AddEnergyBin(TMEBin* e);
    void Add3DBin(TMEBin* e, Int_t nCTbin, Int_t iCTbin, Int_t nMMbin, Int_t iMMbin);
    void AddDataObject(TObject* o);
    void AddEmptyTargetObject(TObject* o);
    void AddSignalObject(TObject* o);
    void AddSignal1Object(TObject* o);
    void AddSignal2Object(TObject* o);
    void AddBackgroundObject(TObject* o);
    void AddTotalObject(TObject* o);
    void AddFitObject(TObject* o);
    void AddArrayObject(Int_t n, TObject** o);

    void GetHistoRange();
    void Draw(const Char_t* szType="all", Double_t xMin = 0., Double_t xMax = 300.);
    void DrawFitCollection(Double_t xMin = -1000., Double_t xMax = 1000.);
    void DrawPaperCollection(Double_t xMin = -1000., Double_t xMax = 1000.);
    void DrawArrayCollection(Int_t iRebin = 1, Double_t xMin = -1000., Double_t xMax = 1000.);

    ClassDef(TMObjectCollection, 1)  // Bin of differential cross section
};

#endif

