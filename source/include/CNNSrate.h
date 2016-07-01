#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

double SMrate(double ErKeV, paramList *pList, int detj);

double intCNNSrate(double Er_min, double Er_max, paramList *pList, int detj);

void SMrateInit( paramList *pList, int detj);

double intSMrate(double Er_min, double Er_max, paramList *pList, int detj);
double diffSMrate(double Er, paramList *pList, int detj);
