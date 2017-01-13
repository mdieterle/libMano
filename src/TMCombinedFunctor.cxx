/*************************************************************************
 * Author: Manuel Dieterle, 2012
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMCombinedFunctor                                                                //
//                                                                      //
// Class for collecting any kind of TObjects                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TMCombinedFunctor.h"

ClassImp(TMCombinedFunctor)

//______________________________________________________________________________
TMCombinedFunctor::TMCombinedFunctor(TH1* hA, TH1* hB, const Char_t* szName, const Char_t* szTitle)
    : TNamed(szName, szTitle)
{
   fA = hA;
   fB = hB;
}

//______________________________________________________________________________
TMCombinedFunctor::~TMCombinedFunctor()
{
    // Destructor.
}

//______________________________________________________________________________
