#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

double sterileRate(double ErKeV, paramList *pList, int detj, int fluxj);

//rates computed with interpolation
double intSterileRate(double Er_min, double Er_max, paramList *pList, int detj);
double diffSterileRate(double Er, paramList *pList, int detj);
