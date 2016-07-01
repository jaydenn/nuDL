#include <iostream>
#include <string.h>
#include <cmath>
#include "DEIntegrator.h"
#include "detectorFunctions.h"
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

//returns integrated total # background events per tonne/year for bg type, with recoil  Er_min < Er < Er_max
double intBgRate(detector det, double Er_min, double Er_max)  						  
{   
    return 0; //gsl_spline_eval_integ(det.background, Er_min, Er_max, det.accel);
}

double diffBgRate(detector det, double Er)
{
    return 0; //gsl_spline_eval(det.background, Er, det.accel);
}

//only need to calculate background once, store in a table for interpolation, stored as N/t/year/keV
int InitializeBackground(detector *det)
{
    //get values of bg at relevant energies
    double Er[500];
    double Bg[500];
    for(int i=0; i<500; i++)
    {
        Er[i] = det->ErL + (double)i*(det->ErU-det->ErL)/499;
        Bg[i] = 0; //detBackground(Er[i], det->bg) * 1000.0 * 365.24;
    }
    
    //create gsl interpolation object
    return gsl_spline_init(det->background,Er,Bg,500);

}

