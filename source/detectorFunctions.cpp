#include <iostream>
#include <cstdio>
#include <cmath>
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif
#ifndef DETECTORSTRUCT_H
	#include "detectorStruct.h"
#endif	

double detEff(double Er, int type)
{
    switch( type ) 
    {
        case 0: 
            return 1;
        case 1: 
            return .5;
        default:
            printf("invalid detector efficiency\n"); 
            return NAN;
    }
}

//return background in events/kg/day/keV
double detBackground(double Er, int type)
{
    switch( type ) 
    {
        case 0: 
            return 1e-99;
        case 1: 
            return 100.0;
	    case 2: 
            return 10.0;
        default:
            printf("invalid detector background\n"); 
            return NAN; 
    }
}

double detRes(double Er, int type)
{
    switch( type ) 
    {
        case 0: 
            return 0;
        case 1: 
            return 0;
        case 2: 
            return 0;
        case 3: 
            return 0;
        case 4: 
            return 0;
        default:
            printf("invalid detector resolution\n"); 
            return NAN; 
    }
}

//returns integrated total # background events per tonne/year for bg type, with recoil  Er_min < Er < Er_max
double intBgRate(detector det, double Er_min, double Er_max)  						  
{   
    return gsl_spline_eval_integ(det.background, Er_min, Er_max, det.accel);
}

double diffBgRate(detector det, double Er)
{
    return gsl_spline_eval(det.background, Er, det.accel);
}
