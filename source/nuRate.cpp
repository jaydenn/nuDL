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

#include "physicalConstants.h";

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
		        - intInvEnu	  * 2*ErGeV*( pow(pList->qA,2)-2*pList->qA*pList->qV + pow(pList->qV,2) ) 
		        + intInvEnuSq * ( ErGeV*ErGeV*( pow(pList->qA,2)-2*pList->qA*pList->qV + pow(pList->qV,2) ) + ErGeV*Mt*(pow(pList->qA,2) - pow(pList->qV,2)) )
	          ); 
    
}

//returns nu scattering rate per target/day/keV with sterile oscillations
double nuRateOsc(double ErKeV, paramList *pList, double Mt, int fluxj)						  
{

	double ErGeV = ErKeV/GeVtoKeV;

    double intConst    = fluxIntegralOsc( ErGeV, pList, Mt,  0, fluxj);
	double intInvEnu   = fluxIntegralOsc( ErGeV, pList, Mt, -1, fluxj);
	double intInvEnuSq = fluxIntegralOsc( ErGeV, pList, Mt, -2, fluxj);

    return pow(GFERMI,2) / ( 2 * M_PI ) * Mt / GeVtoKeV * secsPerDay
	        * ( 
		          intConst	  * 2*( pow(pList->qA,2) + pow(pList->qV,2) )
		        - intInvEnu	  * 2*ErGeV*pow(pList->qA-pList->qV,2) 
		        + intInvEnuSq * ( ErGeV*ErGeV*pow(pList->qA - pList->qV,2) + ErGeV*Mt*(pow(pList->qA,2) - pow(pList->qV,2)) )
	          ); 
    
}

void rateInit( paramList *pList, int detj, int fluxj, double (*rateFunc)(double, paramList *, int, int), gsl_spline *rateSpline)
{
    double ErkeV[INTERP_POINTS];
    double rate[INTERP_POINTS];

    double linStep,logStep; 
    
    logStep = pow((pList->detectors[detj].ErU*1.01)/(pList->detectors[detj].ErL*0.99),1./INTERP_POINTS);
    linStep = (pList->detectors[detj].ErU*1.01-pList->detectors[detj].ErL*0.99)/INTERP_POINTS;

    //always over and undershoot range so that interpolation is well behaved
    ErkeV[0] = 0.99*pList->detectors[detj].ErL;
    rate[0] = rateFunc( (double)ErkeV[0], pList, detj, fluxj);
    
    for( int i=1; i < INTERP_POINTS; i++ )
    {
        
        if(pList->logBins == 0 )//&& ErkeV[i-1] > 5)
            ErkeV[i] = ErkeV[i-1] + linStep;
        else
            ErkeV[i] = ErkeV[i-1] * logStep;
            
        rate[i] = rateFunc( (double)ErkeV[i], pList, detj, fluxj);	
       //std::cout << i << " " << ErkeV[i] << " " << logStep << std::endl;
    }
    //create gsl interpolation object
    gsl_spline_init(rateSpline,ErkeV,rate,INTERP_POINTS);
}


