#include <iostream>
#include <cmath>
#ifndef GSL_INTERP_H
	#include <gsl/gsl_interp.h>
#endif
#include "nuFlux.h"
#ifndef DETECTORSTRUCT_H
	#include "detectorStruct.h"
#endif
#ifndef PARAMETERSTRUCT_H
	#include "parameterStruct.h"
#endif	

const double GFERMI = 1.1664e-5; //GeV^-2
const double HBAR	= 6.5821e-25; //GeV*s
const double GeVtoKeV = 1e6;
const double secsPerDay = 24*60*60; 
const double MN = 0.9383; //mass of nucleon in GeV
const double ME = 0.000510998; //mass of electron in GeV

//returns nu scattering rate per target/day/keV
double nuRate(double ErKeV, paramList *pList, double Mt)						  
{

	double ErGeV = ErKeV/GeVtoKeV;

    double intConst	   = fluxIntegral( ErGeV, pList, Mt, 0);
	double intInvEnu   = fluxIntegral( ErGeV, pList, Mt, -1);
	double intInvEnuSq = fluxIntegral( ErGeV, pList, Mt, -2);
    
    return pow(GFERMI,2) / ( 2 * M_PI ) * Mt / GeVtoKeV * secsPerDay
	        * ( 
		          intConst	  * 2*( pow(pList->qA,2) + pow(pList->qV,2) )
		        - intInvEnu	  * 2*ErGeV*pow(pList->qA-pList->qV,2) 
		        + intInvEnuSq * ( ErGeV*ErGeV*pow(pList->qA - pList->qV,2) + ErGeV*Mt*(pow(pList->qA,2) - pow(pList->qV,2)) )
	          ); 
    
}

void rateInit( paramList *pList, int detj, double (*rateFunc)(double, paramList *, int), gsl_spline *rateSpline)	
{
    double ErkeV[1000];
    double rate[1000];
    for(int i=0; i<1000; i++)
    {
        ErkeV[i] = pList->detectors[detj].ErL + (double)i*(pList->detectors[detj].ErU-pList->detectors[detj].ErL)/900;
        rate[i] = rateFunc( (double)ErkeV[i], pList, detj);	
    }

    //create gsl interpolation object
    gsl_spline_init(rateSpline,ErkeV,rate,1000);
}


