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

#include "physicalConstants.h"

//returns nu scattering rate per target/day/keV
double nuRate(double ErKeV, paramList *pList, double Mt, int fluxj, double dist, int doOsc)		  
{

	double ErGeV = ErKeV/GeVtoKeV;

    double intConst, intInvEnu, intInvEnuSq;
    if( doOsc == 1)
    {
        intConst    = fluxIntegralOsc( ErGeV, pList, Mt,  0, fluxj, dist);
        intInvEnu   = fluxIntegralOsc( ErGeV, pList, Mt, -1, fluxj, dist);
        intInvEnuSq = fluxIntegralOsc( ErGeV, pList, Mt, -2, fluxj, dist);
    }
    else
    {
        intConst    = fluxIntegral( ErGeV, pList, Mt,  0, fluxj, dist);
	    intInvEnu   = fluxIntegral( ErGeV, pList, Mt, -1, fluxj, dist);
	    intInvEnuSq = fluxIntegral( ErGeV, pList, Mt, -2, fluxj, dist);
    }
    
    //return pow(GFERMI,2) / ( M_PI ) * Mt / GeVtoKeV * secsPerDay * (intConst - intInvEnuSq*Mt*ErGeV/2.)*pow(pList->qV,2); //approx calc for testing, could be used to speed up
    return pow(GFERMI,2) / ( 2 * M_PI ) * Mt / GeVtoKeV * secsPerDay
	        * ( 
		          intConst	  * 2*( pow(pList->qA,2) + pow(pList->qV,2) )
		        - intInvEnu	  * 2*ErGeV*( pow(pList->qA,2)-2*pList->qA*pList->qV + pow(pList->qV,2) ) 
		        + intInvEnuSq * (   ErGeV*ErGeV*( pow(pList->qA,2)-2*pList->qA*pList->qV + pow(pList->qV,2) ) + ErGeV*Mt*(pow(pList->qA,2) - pow(pList->qV,2)) )
	          ); 
    
}

double nuRate(double ErKeV, paramList *pList, double Mt, int fluxj, double dist)	
{
     return nuRate(ErKeV, pList, Mt, fluxj, dist, 0);
}


void rateInit( paramList *pList, int detj, int fluxj, double (*rateFunc)(double, paramList *, int, int), gsl_spline *rateSpline)
{
    double ErkeV[INTERP_POINTS];
    double rate[INTERP_POINTS];

    double linStep,logStep; 


    //always over and undershoot range so that interpolation is well behaved
    logStep = pow(pList->detectors[detj].ErU/(0.99*pList->detectors[detj].ErL),1./(INTERP_POINTS-2));
    linStep = (pList->detectors[detj].ErU-0.99*pList->detectors[detj].ErL)/(INTERP_POINTS-2);

    ErkeV[0] = 0.99*pList->detectors[detj].ErL;
    rate[0] = rateFunc( (double)ErkeV[0], pList, detj, fluxj);
    
    for( int i=1; i < INTERP_POINTS; i++ )
    {
        
        if(pList->logBins == 0 )//&& ErkeV[i-1] > 5)
            ErkeV[i] = ErkeV[i-1] + linStep;
        else
            ErkeV[i] = ErkeV[i-1] * logStep;
            
        rate[i] = rateFunc( ErkeV[i], pList, detj, fluxj);	
        //std::cout << i << " " << ErkeV[i] << " " << logStep << " " <<  rate[i] << std::endl;
    }
    
    gsl_spline_init(rateSpline,ErkeV,rate,INTERP_POINTS);

}


