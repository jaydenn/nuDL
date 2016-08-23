#include <cmath>
#include <iostream>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include "parameterStruct.h"
#include "SMrate.h"
#include "BSMrate.h"
#ifndef DETECTORFUNCTIONS_H
	#include "detectorFunctions.h"
#endif

    
//natural log of Poisson dist: gives more accurate values for small probabilities (because of machine precision)
double logPoisson(double obs, double expect)
{
    if ( expect > 0. && obs > 0. )
    {
        if(obs>200)
            return -expect + (double) obs * log( expect ) - ( obs * log( obs ) - obs );  // Sterling's approx.
        else           
            return -expect + (double) obs * log( expect ) - gsl_sf_lngamma( obs+1 );     // or return gsl_ran_poisson_pdf (obs,expect); 
     }
    else
        return -1E299;
}

//likelihood function for binned data
double logLikelihood(paramList *pList)
{
        
	//Calculate log-likelihood
	double loglike = 0;
    double Er_min, Er_max;
    double l,SM,BG,BSM;
    
    //loop over detectors
    for(int detj=0; detj < pList->ndet; detj++)
    {   
        Er_min = pList->detectors[detj].ErL;
        //loop over recoil energy bins
        for(int i=0; i< pList->detectors[detj].nbins; i++)			
        {
            //set bin limits
            Er_max = Er_min + pList->detectors[detj].binW[i];
            
            SM  = intSMrate( Er_min, Er_max, pList, detj);
            BG  = pList->detectors[detj].BgNorm * intBgRate( pList->detectors[detj], Er_min, Er_max);
            BSM = intBSMrate( Er_min, Er_max, pList, detj, pList->signalNorm); 

            l = logPoisson( pList->detectors[detj].binnedData[i], pList->detectors[detj].exposure*(SM+BG+BSM)+1e-99);
            loglike += l;
            //std::cout << " " << i << ": bg " << BG <<  " sm " << SM << " bsm " << BSM << " obs " << pList->detectors[detj].binnedData[i] << " l " << l << " tot " << loglike << std::endl;
            Er_min = Er_max; //update lower bin limit
        } 
        
    }
    
    if( isinf(loglike) )
        return -1e200;
    else
        return loglike;
}

