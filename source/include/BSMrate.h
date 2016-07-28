#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

double BSMrate(double ErKeV, paramList *pList, int detj);

//rates computed with interpolation
double intBSMrate(double Er_min, double Er_max, paramList *pList, int detj, double sigNorm);
double diffBSMrate(double Er, paramList *pList, int detj, double sigNorm);
