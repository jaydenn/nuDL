#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

int nuFluxInit(paramList *pList, std::string sourceName);

double nuFlux(double Enu, paramList *pList);

double fluxIntegral(double ErGeV,  paramList *pList, double Mt, int EnuPow);
