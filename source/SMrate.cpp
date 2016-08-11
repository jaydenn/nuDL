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

//returns SM rate per kg/year/keV
double SMrate(double ErKeV, paramList *pList, int detj)							  
{
    
	double rate = 0;
	paramList pListSM = *pList;
    double targetsPerKG;

    for(int i=0;i<pList->detectors[detj].nIso;i++)
	{
        targetsPerKG = GeVperKG/(MN*pList->detectors[detj].isoA[i]); //how many targets per kg of detector
		
    	//if(pList->nucScat)
	    {
	        pListSM.qA = pList->detectors[detj].isoSN[i]*(-0.427*-0.501163+0.842*0.506875) + pList->detectors[detj].isoSZ[i]*(-0.427*0.506875+0.842*-0.501163);	 
		    pListSM.qV = (- 0.512213 * (pList->detectors[detj].isoA[i] - pList->detectors[detj].isoZ[i]) + .0304*pList->detectors[detj].isoZ[i] )* ffactorSI( pList->detectors[detj].isoA[i], ErKeV);	 
		    rate += targetsPerKG * pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListSM, MN*pList->detectors[detj].isoA[i]);
	    }	
	    //if(pList->elecScat)
	    {
		    pListSM.qA = 0.5;
		    pListSM.qV = 0.5+2*0.2312;
		    int Ne=pList->detectors[detj].isoZ[i];
	        while(pList->detectors[detj].ionization[-1+Ne--] > ErKeV && Ne>0);
		    rate += (Ne * targetsPerKG) * pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListSM, ME);
	    }
		
    }
    
	return rate; 
}

double diffSMrate(double ErkeV, paramList *pList, int detj)						  
{   
    double rate=1e-99;
    for(int i=0; i< pList->source.numFlux; i++)
    {
        if( ErkeV < (pList->detectors[detj].ErL + (double)999*(pList->detectors[detj].ErU-pList->detectors[detj].ErL)/900) )
            rate += gsl_spline_eval(pList->detectors[detj].signalSM[i], ErkeV, pList->detectors[detj].accelSM[i]);
    }
    return rate;
}

double intSMrate(double Er_min, double Er_max, paramList *pList, int detj)						  
{   
    double rate = 0;
    for(int i=0; i< pList->source.numFlux; i++)
        rate += gsl_spline_eval_integ(pList->detectors[detj].signalSM[i], Er_min, Er_max, pList->detectors[detj].accelSM[i]);
    
    return rate;
}

