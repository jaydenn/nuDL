#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

int nuFluxInit(paramList *pList, std::string sourceName);

double nuFlux(double Enu, paramList *pList, int fluxj);

double fluxIntegral(double ErGeV,  paramList *pList, double Mt, int EnuPow, int fluxj, double dist);

double fluxIntegralOsc(double ErGeV,  paramList *pList, double Mt, int EnuPow, int fluxj, double dist);
