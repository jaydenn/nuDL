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

const double MN = 0.9383; //mass of nucleon in GeV
const double ME = 0.000510998; //mass of electron in GeV
const double GeVperKG = 5.6094e26;
const double SSW = 0.2387; //sin^2(theta_w)

//returns rate per kg/year/keV for the jth flux
double sterileRate(double ErKeV, paramList *pList, int detj, int fluxj)							  
{
    
	double rate = 0;
	paramList pListSM = *pList;
    double targetsPerKG;

    for(int i=0;i<pList->detectors[detj].nIso;i++)
	{
        targetsPerKG = GeVperKG/(MN*pList->detectors[detj].isoA[i]); //how many targets per kg of detector
		
    	if(pList->nucScat)
	    {
	        pListSM.qA = pList->detectors[detj].isoSN[i]*(-0.427*-0.501163+0.842*0.506875) + pList->detectors[detj].isoSZ[i]*(-0.427*0.506875+0.842*-0.501163);	 
		    pListSM.qV = (- 0.512213 * (pList->detectors[detj].isoA[i] - pList->detectors[detj].isoZ[i]) + (1-4*SSW)*pList->detectors[detj].isoZ[i] )* ffactorSI( pList->detectors[detj].isoA[i], ErKeV);	 
		    rate += targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListSM, MN*pList->detectors[detj].isoA[i], fluxj);
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
		        rate += pList->source.survProb[fluxj] * (pList->detectors[detj].isoZ[i]-Ne) * targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRateOsc( ErKeV, &pListSM, ME, fluxj);

		        pListSM.qA = -0.5;
		        pListSM.qV = 2*SSW-0.5;
		        rate += (1-pList->source.survProb[fluxj]) * (pList->detectors[detj].isoZ[i]-Ne) * targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRateOsc( ErKeV, &pListSM, ME, fluxj);
		    }
		    else
		    {
		        pListSM.qA = -0.5;
		        pListSM.qV = 0.5+2*SSW;
		        rate += (pList->detectors[detj].isoZ[i]-Ne) * targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRateOsc( ErKeV, &pListSM, ME, fluxj);
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
        if( ErkeV < (pList->detectors[detj].ErL + (double)999*(pList->detectors[detj].ErU-pList->detectors[detj].ErL)/900) )
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

