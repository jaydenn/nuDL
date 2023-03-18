#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

double NSIrate(double ErKeV, paramList *pList, int detj, int fluxj);

//rates computed with interpolation
double intNSIrate(double Er_min, double Er_max, paramList *pList, int detj);
double diffNSIrate(double Er, paramList *pList, int detj);
