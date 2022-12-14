#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

void rateInit(paramList *pList, int detj, int fluxj, double (*rateFunc)(double, paramList *, int, int), gsl_spline *rateSpline);

double nuRate(double ErKeV, paramList *pList, double Mt, int fluxj, double dist);
double nuRate(double ErKeV, paramList *pList, double Mt, int fluxj, double dist, int doOsc);
