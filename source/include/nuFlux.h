#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

int initFlux(paramList *pList);

double nuFlux(double Enu, paramList *pList);

double fluxIntegral(double ErGeV,  paramList *pList, double Mt, int EnuPow);
