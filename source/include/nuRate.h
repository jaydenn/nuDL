#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

void rateInit(paramList *pList, int detj, double (*rateFunc)(double, paramList *, int), gsl_spline *rateSpline);

double nuRate(double ErKeV, paramList *pList, double Mt);
