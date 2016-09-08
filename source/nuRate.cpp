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

    for(int i=0; i<5000; i++)
    {
        //always over and undershoot range so that interpolation is well behaved
        if(pList->logBins == 1)
            ErkeV[i] = pow(10, log10(0.95*pList->detectors[detj].ErL) + (double)i*(log10(pList->detectors[detj].ErU)-log10(pList->detectors[detj].ErL))/4800);
        else
            ErkeV[i] = 0.95*pList->detectors[detj].ErL + (double)i*(pList->detectors[detj].ErU-pList->detectors[detj].ErL)/4800;
            
        rate[i] = rateFunc( (double)ErkeV[i], pList, detj, fluxj);	
    }
    //create gsl interpolation object
    gsl_spline_init(rateSpline,ErkeV,rate,5000);
}


