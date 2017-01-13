
/*************************************************************************
 * Author: Manuel Dieterle, 2012
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMCollection                                                         //
//                                                                      //
// Class for a collection of TObjects                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef myTMCollection
#define myTMCollection

#include "TList.h"
#include "TNamed.h"
#include "TF1.h"
#include "TH1.h"
#include "TMEBin.h"
#include "TCanvas.h"

class TMCollection : public TNamed
{

protected:
    TList* fCollection;				// list of objects in collection

public:

//    TMCollection() : TNamed(), fCollection(0) { }
    TMCollection(const Char_t* szName = "", const Char_t* szTitle = "");
    virtual ~TMCollection();

    TList* GetCollection() const { return fCollection; }

    void AddObject(TObject* o);
    void ListObjects();
    TObject* GetObject(const Char_t* szName);
    Bool_t IsObjectType(const Char_t* szName, const Char_t* szType);
    TF1* GetTF1(const Char_t* szName);
    TH1* GetTH1(const Char_t* szName);

    ClassDef(TMCollection, 1)  // Bin of differential cross section
};

#endif

