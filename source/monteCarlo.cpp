#include <cmath>
#include <iostream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <iomanip>
#include "parameterStruct.h"
#ifndef DETECTORFUNCTIONS_H
	#include "detectorFunctions.h"
#endif
int SEED=0;

/*
//Generates random number of randomly distributed according to dN/dE for parameters in cube, for detector det (0 or 1, corresponding to order in sampling.par), recoil energies [keV] are stored in MCdata
int generateUnbinnedData(WIMPpars W, detector *det, int b, int simSeed)
{

    double signal = intWIMPrate( det->ErL, det->ErU, W, *det) * det->exposure;   
    double background= b*intBgRate(det, det->ErL, det->ErU) * det->exposure; 
     
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T=gsl_rng_default;
    r=gsl_rng_alloc(T);
    gsl_rng_set(r, simSeed + SEED++);
    
    det->nEvents = gsl_ran_poisson(r,signal+background);
    
    det->unbinnedData = new double[(int)det->nEvents];
    
    double norm = diffWIMPrate( det->ErL, W, *det) + b*diffBgRate(*det,det->ErL);
    int i = 0;
    double Er,y;
    while ( i < det->nEvents)
    {
        Er = gsl_rng_uniform(r)*(det->ErU-det->ErL) + det->ErL;
        y = gsl_rng_uniform(r);
        if( (diffWIMPrate( det->ErL, W, *det) + b*diffBgRate(*det,det->ErL))/norm  > y )
        {
            det->unbinnedData[i] = Er;
            i++;
        }
    }
    
    return 0;
}*/

int generateBinnedData(int detj, int b, int simSeed)
{

    double Er_min, Er_max;
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T=gsl_rng_default;
    r=gsl_rng_alloc(T);
    gsl_rng_set(r, simSeed + SEED++);

    //total signal and background, for setting bin size
double signal = 0;//intWIMPrate( det->ErL, det->ErU, W, det) * det->exposure;
double background = 0;//b * intBgRate(det, det->ErL, det->ErU) * det->exposure; 
    //cout << signal << " " << background << endl;
    //somewhat arbitrary choice of number of bins.. seems to work for exponential data
    
    //setup bins
    /*
    det->nbins = floor( sqrt(signal+background) )+2;
    det->binW = ( det->ErU - det->ErL ) / ( (double) det->nbins);

    try
    {
        //det->binnedData = new double[det->nbins];
    }
    catch (std::bad_alloc& ba)
    {
      std::cerr << "bad_alloc caught: " << ba.what() << std::endl << "you requested: " << det->nbins << " doubles" <<std::endl;
      return 1;
    }
    */
    
/*    for(int i=0; i<det->nbins; i++)
    {
        Er_min = (double)i*det->binW+det->ErL;
        Er_max = (double)(i+1)*det->binW+det->ErL;
         
//        background = b * intBgRate(det, Er_min, Er_max) * det->exposure;
//        signal = intWIMPrate( Er_min, Er_max, W, det) * det->exposure; 

        if( W->asimov == 1) 
            det->binnedData[i] = gsl_ran_poisson(r,signal+background);
        else
            det->binnedData[i] = signal + background;            
        
        det->nEvents += det->binnedData[i];
    }
  */  
    return 0;
     
}

