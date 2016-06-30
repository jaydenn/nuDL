#include <cmath>
#include <iostream>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include "parameterStruct.h"
#include "CNNSrate.h"
#include "detectors.h"
#include "detectorFunctions.h"
#include <assert.h>
    
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
	double signal, background, l;
	double loglike = 0;
    double Er_min, Er_max;

    //loop over detectors
    for(int j=0; j < pList->ndet; j++)
    {   
        //loop over recoil energy bins
        for(int i=0; i< pList->detectors[j].nbins; i++)			
        {
            //set bin limits
            Er_min = (double)i*pList->detectors[j].binW + pList->detectors[j].ErL;
            Er_max = (double)(i+1)*pList->detectors[j].binW + pList->detectors[j].ErL;
            
           // background = intBgRate(  Er_min, Er_max, pList, j) * pList->detectors[j].exposure;
           // signal     = intCNNSrate(Er_min, Er_max, pList, j) * pList->detectors[j].exposure; 
            
            l = logPoisson( pList->detectors[j].binnedData[i], signal+background+1e-99);
            loglike += l;
            //cout << " " << i << ": bg " << background <<  " s " << signal << " obs " << pL->detectors[j].binnedData[i] << " l " << l << " tot " << loglike << endl;
        } 
        
    }
   
    if( isinf(loglike) )
        return -1e200;
    else
        return loglike;
}

