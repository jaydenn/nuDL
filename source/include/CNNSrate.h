#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

double SMrate(double ErKeV, paramList *pList, int detj);
double BSMrate(double ErKeV, paramList *pList, int detj);

double intCNNSrate(double Er_min, double Er_max, paramList *pList, int detj);

void rateInit(paramList *pList, int detj, double (*rateFunc)(double, paramList *, int), gsl_spline *rateSpline);

double intSMrate(double Er_min, double Er_max, paramList *pList, int detj);
double diffSMrate(double Er, paramList *pList, int detj);
double intBSMrate(double Er_min, double Er_max, paramList *pList, int detj);
double diffBSMrate(double Er, paramList *pList, int detj);
