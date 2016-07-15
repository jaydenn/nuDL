#include <cmath>
#include <iostream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <iomanip>
#include "SMrate.h"
#include "BSMrate.h"
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
    
    //total standard model, beyond SM and background rates
    double SM, BSM, BG;
    if(pList->BSM)
        BSM = pList->nuFluxNorm * pList->signalNorm * intBSMrate( pList->detectors[detj].ErL, pList->detectors[detj].ErU, pList, detj);
    else
        BSM = 0;
        
    SM = pList->nuFluxNorm * intSMrate( pList->detectors[detj].ErL, pList->detectors[detj].ErU, pList, detj);
    BG = pList->detectors[detj].BgNorm * b * intBgRate(pList->detectors[detj], pList->detectors[detj].ErL, pList->detectors[detj].ErU) ;     

    //std::cout << SM << " " << BG << " " << BSM << std::endl;

    //setup bins ~somewhat arbitrary choice of number of bins.. seems to work for exponential data
    pList->detectors[detj].nbins = floor( sqrt(  pList->detectors[detj].exposure * ( SM + BSM + BG ) ) ) + 2;

    if (pList->detectors[detj].nbins > 50)
        pList->detectors[detj].nbins = 50;
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
        
        if(pList->BSM)
            BSM = pList->nuFluxNorm * pList->signalNorm * intBSMrate( Er_min, Er_max, pList, detj);
        else
            BSM = 0;
            
        SM = pList->nuFluxNorm * intSMrate( Er_min, Er_max, pList, detj); 
        BG = pList->detectors[detj].BgNorm * b * intBgRate(pList->detectors[detj], Er_min, Er_max) ;
        
        if( pList->asimov == 1) 
            pList->detectors[detj].binnedData[i] = pList->detectors[detj].exposure *( SM + BSM + BG );
        else
            pList->detectors[detj].binnedData[i] = gsl_ran_poisson(r, pList->detectors[detj].exposure *( SM + BSM + BG ));                       
        
        pList->detectors[detj].nEvents += pList->detectors[detj].binnedData[i];
    }
   
    return 0;
     
}

