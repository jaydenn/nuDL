#include <iostream>
#include <cmath>
#ifndef GSL_INTERP_H
	#include <gsl/gsl_interp.h>
#endif
#include "formFactorSI.h"
#include "nuRate.h"
#ifndef DETECTORSTRUCT_H
	#include "detectorStruct.h"
#endif
#ifndef PARAMETERSTRUCT_H
	#include "parameterStruct.h"
#endif	

#include "physicalConstants.h"

//returns rate with SM couplings but including sterile oscillations per kg/year/keV for the jth flux
double sterileRate(double ErKeV, paramList *pList, int detj, int fluxj)
{
    
	double rate = 0;
	paramList pListSM = *pList;
    double targetsPerKG;

    for(int i=0;i<pList->detectors[detj].nIso;i++)
	{
        targetsPerKG =  GeVperKG/(AMU*pList->detectors[detj].AM); //how many targets per kg of detector
		
    	if(pList->nucScat)
	    {
	        pListSM.qA = 4.0/3.0 * (pList->detectors[detj].isoJN[i]+1) / pList->detectors[detj].isoJN[i] * ( pList->detectors[detj].isoSN[i]*GAN + pList->detectors[detj].isoSZ[i]*GAP );	 
		    pListSM.qV = ( GVN * (pList->detectors[detj].isoA[i] - pList->detectors[detj].isoZ[i]) + GVP * pList->detectors[detj].isoZ[i] )* ffactorSIhelm( pList->detectors[detj].isoA[i], sqrt(2*pList->detectors[detj].isoA[i]*AMU*ErKeV/GeVtoKeV) );	
		    rate += targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListSM, AMU*pList->detectors[detj].isoA[i], fluxj, 1);
	    }	
	    if(pList->elecScat)
	    {
	        int Ne=0;
	        while(pList->detectors[detj].ionization[i][Ne] > ErKeV && Ne < pList->detectors[detj].isoZ[i]) 
	            Ne++;
	    
	        if(pList->source.isSolar[fluxj] == 1)
	        {
		        pListSM.qA = 0.5;
		        pListSM.qV = 2*SSW+0.5;
		        rate += pList->source.survProb[fluxj] * (pList->detectors[detj].isoZ[i]-Ne) * targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListSM, ME, fluxj, 1);

		        pListSM.qA = -0.5;
		        pListSM.qV = 2*SSW-0.5;
		        rate += (1-pList->source.survProb[fluxj]) * (pList->detectors[detj].isoZ[i]-Ne) * targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListSM, ME, fluxj, 1);
		    }
		    else
		    {
		        pListSM.qA = -0.5;
		        pListSM.qV = 0.5+2*SSW;
		        rate += (pList->detectors[detj].isoZ[i]-Ne) * targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListSM, ME, fluxj, 1);
		    }
		    
	    }
		
    }
    
	return rate; 
}

double diffSterileRate(double ErkeV, paramList *pList, int detj)						  
{   
    double rate=1e-99;
    for(int i=0; i< pList->source.numFlux; i++)
    {
        if( ErkeV >= pList->detectors[detj].ErL || ErkeV <= pList->detectors[detj].ErU )
            rate += pList->source.nuFluxNorm[i] * gsl_spline_eval(pList->detectors[detj].signalBSM1[i], ErkeV, pList->detectors[detj].accelBSM1[i]);
    }
    return rate;
}

double intSterileRate(double Er_min, double Er_max, paramList *pList, int detj)						  
{   
    double rate = 0;
    for(int i=0; i< pList->source.numFlux; i++)
        rate += pList->source.nuFluxNorm[i] * gsl_spline_eval_integ(pList->detectors[detj].signalBSM1[i], Er_min, Er_max, pList->detectors[detj].accelBSM1[i]);
    
    return rate;
}

