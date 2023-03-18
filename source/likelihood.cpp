#include <cmath>
#include <iostream>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include "parameterStruct.h"
#include "nuRate.h"
#include "SMrate.h"
#include "BSMrate.h"
#include "sterileRate.h"
#ifndef DETECTORFUNCTIONS_H
	#include "detectorFunctions.h"
#endif

    
//natural log of Poisson dist: gives more accurate values for small probabilities (because of machine precision)
double logPoisson(double obs, double expect)
{
    if ( expect > 0. && obs >= 0. )
        return -expect + obs * log( expect ) - gsl_sf_lngamma( obs+1 );
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

            l = logPoisson( pList->detectors[detj].binnedData[i], pList->detectors[detj].exposure*(SM+BG+BSM)+1e-200);
            loglike += l;
            Er_min = Er_max; //update lower bin limit
            //std::cout << "l: " << i << " " << SM << " " << BG << " " << BSM << " " << pList->detectors[detj].binnedData[i] << " " << l << std::endl;
        } 
        
    }
    if( std::isinf(loglike) )
        return -1e200;
    else
        return loglike;
}

//likelihood function for binned data
double logLikelihoodSM(paramList *pList)
{
        
	//Calculate log-likelihood
	double loglike = 0;
    double Er_min, Er_max;
    double l,SM,BG;
    
    //loop over detectors
    for(int detj=0; detj < pList->ndet; detj++)
    {   
        Er_min = pList->detectors[detj].ErL;
        //loop over recoil energy bins
        for(int i=0; i< pList->detectors[detj].nbins; i++)			
        {
            //set bin limits
            Er_max = Er_min + pList->detectors[detj].binW[i];
            SM  = pList->signalNorm * intSMrate( Er_min, Er_max, pList, detj);
            BG  = pList->detectors[detj].BgNorm * intBgRate( pList->detectors[detj], Er_min, Er_max);
            
            l = logPoisson( pList->detectors[detj].binnedData[i], pList->detectors[detj].exposure*(SM+BG));
            loglike += l;
            //std::cout << " SM" << i << ": bg " << pList->detectors[detj].exposure*BG <<  " sm " << pList->detectors[detj].exposure*SM <<  " exp " << pList->detectors[detj].exposure*(BG+SM) << " obs " << pList->detectors[detj].binnedData[i] << " l " << l << " tot " << loglike << std::endl;
            Er_min = Er_max; //update lower bin limit
        } 
    }

    if( std::isinf(loglike) )
        return -1e200;
    else
        return loglike;
}

double logLikelihoodSterile(paramList *pList)
{

	//Calculate log-likelihood
	double loglike = 0;
    double Er_min, Er_max;
    double l,SM,BG;
    
    //loop over detectors
    for(int detj=0; detj < pList->ndet; detj++)
    {   
        Er_min = pList->detectors[detj].ErL;
        //loop over recoil energy bins
        for(int i=0; i<pList->detectors[detj].nbins; i++)			
        {
            //set bin limits
            Er_max = Er_min + pList->detectors[detj].binW[i];
            SM  = intSterileRate( Er_min, Er_max, pList, detj);
            BG  = pList->detectors[detj].BgNorm * intBgRate( pList->detectors[detj], Er_min, Er_max);

            l = logPoisson( pList->detectors[detj].binnedData[i], pList->detectors[detj].exposure*(SM+BG));
            loglike += l;
            //std::cout << " " << i << ": bg " << pList->detectors[detj].exposure*BG <<  " sm " << pList->detectors[detj].exposure*SM <<  " exp " << pList->detectors[detj].exposure*(BG+SM) << " obs " << pList->detectors[detj].binnedData[i] << " l " << l << " tot " << loglike << std::endl;
            Er_min = Er_max; //update lower bin limit
        } 
    }

    if( std::isinf(loglike) || pList->ss2Theta14 > 1 || pList->ss2Theta14 < 0 )
        return -1e200;
    else
        return loglike;
}

double logLikelihoodSterile(paramList *pList, double ss2t, double delMsq)
{

	//Calculate log-likelihood
	double loglike = 0;
    double Er_min, Er_max;
    double l,sig,BG;
    
    pList->delMsqGeV = delMsq;
    pList->ss2Theta14 = ss2t;
    
    //loop over detectors
    for(int detj=0; detj < pList->ndet; detj++)
    {   

        for( int fluxj=0; fluxj < pList->source.numFlux; fluxj++)
            rateInit( pList, detj, fluxj, &sterileRate, pList->detectors[detj].signalBSM1[0]);
        
        Er_min = pList->detectors[detj].ErL;
        //loop over recoil energy bins
        for(int i=0; i<pList->detectors[detj].nbins; i++)			
        {
            //set bin limits
            Er_max = Er_min + pList->detectors[detj].binW[i];
            sig  = intSterileRate( Er_min, Er_max, pList, detj);
            BG  = pList->detectors[detj].BgNorm * intBgRate( pList->detectors[detj], Er_min, Er_max);

            l = logPoisson( pList->detectors[detj].binnedData[i], pList->detectors[detj].exposure*(sig+BG));
            loglike += l;
            //std::cout << " " << i << ": bg " << pList->detectors[detj].exposure*BG <<  " sm " << pList->detectors[detj].exposure*SM <<  " exp " << pList->detectors[detj].exposure*(BG+SM) << " obs " << pList->detectors[detj].binnedData[i] << " l " << l << " tot " << loglike << std::endl;
            Er_min = Er_max; //update lower bin limit
        } 
    }

    if( std::isinf(loglike) || pList->ss2Theta14 > 1 || pList->ss2Theta14 < 0 )
        return -1e200;
    else
        return loglike;
}