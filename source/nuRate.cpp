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
double nuRate(double ErKeV, paramList *pList, double Mt, int fluxj)						  
{

	double ErGeV = ErKeV/GeVtoKeV;

    double intConst    = fluxIntegral( ErGeV, pList, Mt,  0, fluxj);
	double intInvEnu   = fluxIntegral( ErGeV, pList, Mt, -1, fluxj);
	double intInvEnuSq = fluxIntegral( ErGeV, pList, Mt, -2, fluxj);
    
    return pow(GFERMI,2) / ( 2 * M_PI ) * Mt / GeVtoKeV * secsPerDay
	        * ( 
		          intConst	  * 2*( pow(pList->qA,2) + pow(pList->qV,2) )
		        - intInvEnu	  * 2*ErGeV*pow(pList->qA-pList->qV,2) 
		        + intInvEnuSq * ( ErGeV*ErGeV*pow(pList->qA - pList->qV,2) + ErGeV*Mt*(pow(pList->qA,2) - pow(pList->qV,2)) )
	          ); 
    
}

void rateInit( paramList *pList, int detj, int fluxj, double (*rateFunc)(double, paramList *, int, int), gsl_spline *rateSpline)
{
    double ErkeV[5000];
    double rate[5000];

    double linStep,logStep; 
    
    if(pList->logBins == 0)
    {
        logStep = pow(pList->detectors[detj].ErU/pList->detectors[detj].ErL,1/4800.0);
        linStep = (pList->detectors[detj].ErU-5)/(4800 - log(5/pList->detectors[detj].ErL)/log(logStep) );
    }
    else
        logStep = pow(pList->detectors[detj].ErU/pList->detectors[detj].ErL,1/4800.0);
        linStep = (pList->detectors[detj].ErU-5)/(4800 - log(5/pList->detectors[detj].ErL)/log(logStep) );
    
    ErkeV[0] = 0.95*pList->detectors[detj].ErL; 
    rate[0] = rateFunc( (double)ErkeV[0], pList, detj, fluxj);
    
    for( int i=1; i < 5000; i++ )
    {
        //always over and undershoot range so that interpolation is well behaved
        if(pList->logBins == 0 && ErkeV[i-1] > 5)
            ErkeV[i] = ErkeV[i-1] + linStep;
        else
            ErkeV[i] = ErkeV[i-1] * logStep;
            
        rate[i] = rateFunc( (double)ErkeV[i], pList, detj, fluxj);	
      // std::cout << i << " " << ErkeV[i] << std::endl;
    }
    //create gsl interpolation object
    gsl_spline_init(rateSpline,ErkeV,rate,5000);
}


