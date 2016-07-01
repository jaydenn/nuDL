#include <cmath>
#include <iostream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <iomanip>
#include "CNNSrate.h"
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif
#ifndef DETECTORFUNCTIONS_H
	#include "detectorFunctions.h"
#endif

int SEED=0;

int generateBinnedData(paramList *pList, int detj, int b, int simSeed)
{

    double Er_min, Er_max;
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T=gsl_rng_default;
    r=gsl_rng_alloc(T);
    gsl_rng_set(r, simSeed + SEED++);

    //total signal and background, for setting bin size
    double signal     = pList->nuFluxNorm * pList->signalNorm * pList->rateFunc( pList->detectors[detj].ErL, pList->detectors[detj].ErU, pList, detj) * pList->detectors[detj].exposure;
    double background = pList->detectors[detj].BgNorm * b * intBgRate(pList->detectors[detj], pList->detectors[detj].ErL, pList->detectors[detj].ErU) * pList->detectors[detj].exposure;     
    
    //setup bins //somewhat arbitrary choice of number of bins.. seems to work for exponential data
    pList->detectors[detj].nbins = floor( sqrt(signal+background) )+2;
    pList->detectors[detj].binW  = ( pList->detectors[detj].ErU - pList->detectors[detj].ErL ) / ( (double) pList->detectors[detj].nbins);

    try
    {
        pList->detectors[detj].binnedData = new double[pList->detectors[detj].nbins];
    }
    catch (std::bad_alloc& ba)
    {
      std::cerr << "bad_alloc caught: " << ba.what() << std::endl << "you requested: " << pList->detectors[detj].nbins << " doubles" <<std::endl;
      return 1;
    }

    
    for(int i=0; i<pList->detectors[detj].nbins; i++)
    {
        Er_min = (double)i*pList->detectors[detj].binW+pList->detectors[detj].ErL;
        Er_max = (double)(i+1)*pList->detectors[detj].binW+pList->detectors[detj].ErL;
        
        background = pList->detectors[detj].BgNorm * b * intBgRate(pList->detectors[detj], Er_min, Er_max) * pList->detectors[detj].exposure;
        signal     = pList->nuFluxNorm * pList->signalNorm * pList->rateFunc( Er_min, Er_max, pList, detj) * pList->detectors[detj].exposure; 

        if( pList->asimov == 1) 
            pList->detectors[detj].binnedData[i] = signal + background;
        else
            pList->detectors[detj].binnedData[i] = gsl_ran_poisson(r,signal+background);                       
        
        pList->detectors[detj].nEvents += pList->detectors[detj].binnedData[i];
    }
   
    return 0;
     
}

