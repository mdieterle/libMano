
/*************************************************************************
 * Author: Manuel Dieterle, 2012
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// LinkDef.h                                                            //
//                                                                      //
// libMano dictionary header file.                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifdef __CINT__ 

// turn everything off
#pragma link off all globals; 
#pragma link off all classes; 
#pragma link off all functions; 
#pragma link off all typedef; 

#pragma link C++ nestedclasses; 
#pragma link C++ nestedtypedef; 

#pragma link C++ namespace TMTools;
#pragma link C++ namespace TMCuts;
#pragma link C++ namespace TMFunctor;
#pragma link C++ class TMEBin+; 
#pragma link C++ class TMObjectCollection+; 
#pragma link C++ class TMFitCollection-; 
#pragma link C++ class TMCollection+;
#pragma link C++ class TMFit+;
#pragma link C++ class TMCombinedFit+;
#pragma link C++ class TMCombinedFunctor+;
#pragma link C++ class TMPluto+;
#pragma link C++ class TMEnergyThetaCut+;

#endif

