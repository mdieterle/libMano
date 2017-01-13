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


#include "TMCollection.h"

ClassImp(TMCollection)

//______________________________________________________________________________
TMCollection::TMCollection(const Char_t* szName, const Char_t* szTitle)
    : TNamed(szName, szTitle)
{
   fCollection = new TList();
   fCollection->SetOwner(kTRUE);
}

//______________________________________________________________________________
TMCollection::~TMCollection()
{
   if (fCollection) delete fCollection;
}

//______________________________________________________________________________
void TMCollection::AddObject(TObject* o)
{
   fCollection->Add(o);
   return;
}

//______________________________________________________________________________
void TMCollection::ListObjects()
{
    TIter next(fCollection);
    TObject* o;
    printf("Collection content:\n");
    printf("TYPE            NAME\n");
    while ((o = (TObject*)next()))
       printf("%-15s %s\n", o->ClassName(), o->GetName());

    return;
}

//______________________________________________________________________________
TObject* TMCollection::GetObject(const Char_t* szName)
{
   TObject* obj = fCollection->FindObject(szName);
   return obj ? (TObject*) obj : 0;
}

//______________________________________________________________________________
Bool_t TMCollection::IsObjectType(const Char_t* szName, const Char_t* szType)
{
   if (!strcmp(GetObject(szName)->ClassName(),szType))
      return kTRUE;
   else
      return kFALSE;
}
//______________________________________________________________________________


