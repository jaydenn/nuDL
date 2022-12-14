#ifndef DETECTORSTRUCT_H
	#include "detectorStruct.h"
#endif
#define DETECTORFUNCTIONS_H
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

//detector efficiency as a function of recoil energy [keV]
double detEff(double Er, int type);

//detector background in dru as a function of recoil energy [keV]
double detBackground(double Er, paramList *pList, int detj, int fluxj);		 

//detector background in keV as a function of recoil energy [keV]
double detRes(double Er, int type);

//detector background integrated between two energies in events/t/year
double intBgRate(detector det, double Er_min, double Er_max);		  

//differential background rate /keV/t/year
double diffBgRate(detector det, double Er);							   

//for declaring a new detector
int newDetector(paramList *pList, char *name, double exp, double dist);
