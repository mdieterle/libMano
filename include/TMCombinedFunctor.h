/*************************************************************************
 * Author: Manuel Dieterle, 2012
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMCombinedFunctor                                                            //
//                                                                      //
// Class containing all functions.                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef myTMCombinedFunctor
#define myTMCombinedFunctor

#include "TNamed.h"
#include "TH1.h"

class TMCombinedFunctor : public TNamed
{

protected:

   TH1* fA;
   TH1* fB;

public:
   TMCombinedFunctor() : TNamed(), fA(0), fB(0) {}

   TMCombinedFunctor(TH1* hA, TH1* hB, const Char_t* szName, const Char_t* szTitle);
   virtual ~TMCombinedFunctor();

   Double_t MyFitFunc (Double_t *x, Double_t *par)
   {
      Double_t vA = fA->GetBinContent( fA->GetXaxis()->FindBin( x[0] ) );
      Double_t vB = fB->GetBinContent( fB->GetXaxis()->FindBin( x[0] ) );

      Double_t val = par[0]*vA + par[1]*vB;

      return val;
   };

   ClassDef(TMCombinedFunctor, 0);
};

#endif

