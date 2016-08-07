#include <iostream>
#include <cmath>
#ifndef GSL_INTERP_H
	#include <gsl/gsl_interp.h>
#endif
#include "formFactorSI.h"
#include "nuRate.h"
#include "nuFlux.h"
#ifndef DETECTORSTRUCT_H
	#include "detectorStruct.h"
#endif
#ifndef PARAMETERSTRUCT_H
	#include "parameterStruct.h"
#endif	

const double GFERMI = 1.1664e-5; //GeV^-2
const double HBAR	= 6.5821e-25; //GeV*s
const double GeVtoKeV = 1e6;
const double GeVperKG = 5.6094e26;
const double MN = 0.9383; //mass of nucleon in GeV
const double ME = 0.000510998; //mass of electron in GeV

//returns simplified model light mediator rate per electron/day/keV
double BSMrateE(double ErKeV, paramList *pList, double Mt)						  
{

	double ErGeV = ErKeV/GeVtoKeV;
	double secsPerDay = 24*60*60; 
	
    double convFactor =  1/GeVtoKeV * secsPerDay; //units GeV*s/keV/kg/day
    
    //electron couplings
    double gv = pList->qV;
    double ga = pList->qA;
    
    //makes reading calculations more readable
    double gEs = pList->gEs;   //scalar electron
    double gEp = pList->gEp;   //pseudoscalar electron
    double gNuV = pList->gNuV;  //vector neutrino
    double gNuS = pList->gNuS;  //vector neutrino
    double gEv =  pList->gEv;  //vector electron
    double gEa =  pList->gEa;  //axial-vector electron
    
    double intConst, intInvEnuSq;
        
    intConst	= fluxIntegral( ErGeV, pList, ME, 0);  //units GeV^2/s
    intInvEnuSq = fluxIntegral( ErGeV, pList, ME, -2); //units 1/s

    switch( pList->BSM )
    {
        //scalar
        case 1:
        {
            return convFactor * pow(gEs*gNuS,2) / ( 4 * M_PI ) * pow(ME,2)
                * intInvEnuSq * ErGeV / pow( 2*ErGeV*ME + pow(pList->mMed,2) ,2);
        }
        //pseudoscalar
        case 2:
        {
            return convFactor * pow(gEp*gNuS,2) / ( 8 * M_PI ) * ME
                * intInvEnuSq * pow(ErGeV,2) / pow( 2*ErGeV*ME + pow(pList->mMed,2) ,2);
        }
        //vector
        case 3:
        {
            return convFactor * ( pList->SMinterference1 * sqrt(2)*GFERMI*gv*gNuV*gEv*ME /  M_PI 
                * intConst / ( 2*ErGeV*ME + pow(pList->mMed,2) )
                + pList->SMinterference2 * pow(gEv*gNuV,2) * ME / ( 2 * M_PI ) 
                * intConst / pow( 2*ErGeV*ME + pow(pList->mMed,2) ,2) );   
        }
        //axialvector
        case 4:
        {
            return convFactor * ( pList->SMinterference1 * sqrt(2)*GFERMI*ME*ga*gNuV*gEa /  M_PI 
                * intConst / ( 2*ErGeV*ME + pow(pList->mMed,2) )
                + pList->SMinterference2 * pow(gEa*gNuV,2)* ME / ( 2 * M_PI ) 
                * intConst / pow( 2*ErGeV*ME + pow(pList->mMed,2) ,2) );   
        }
        default:
        {
            return 1e-99;
        }
    }
}

//returns simplified model light mediator rate per nuclei/day/keV
double BSMrateN(double ErKeV, paramList *pList, double Mt)						  
{

	double ErGeV = ErKeV/GeVtoKeV;
	double secsPerDay = 24*60*60; 
	
    double convFactor =  1/GeVtoKeV * secsPerDay; //units GeV*s/keV/kg/day
    
    //nuclear couplings
    double Qv = pList->Qv; //SM vector
    double Qa = pList->Qa; //SM axial coupling
    
    //makes calculations more readable
    double Qs = pList->Qs;  //scalar
    double Qvp = pList->Qvp; //vector
    double Qap = pList->Qap; //axial 

    double intConst, intInvEnu, intInvEnuSq;
    
    intConst	= fluxIntegral( ErGeV, pList, Mt,  0); //units GeV^2/s
    intInvEnu   = fluxIntegral( ErGeV, pList, Mt, -1); //units GeV/s
    intInvEnuSq = fluxIntegral( ErGeV, pList, Mt, -2); //units 1/s

    switch( pList->BSM ) 
    {
        //scalar
        case 1:
        {
            return convFactor * pow(Qs,2) / ( 4 * M_PI ) * pow(Mt,2)
	            * intInvEnuSq * ErGeV / pow( 2*ErGeV*Mt + pow(pList->mMed,2) ,2);
        }
        //pseudoscalar
        case 2:
        {
            return 0;
        }
        //vector
        case 3:
        {
            return convFactor * ( 
                     - pList->SMinterference1 * GFERMI*Mt*Qv*Qvp / (2*sqrt(2)*M_PI * ( 2*ErGeV*Mt + pow(pList->mMed,2) ) ) 
	                     * ( 2*intConst - intInvEnuSq * ErGeV * Mt ) 
	                 + pList->SMinterference2 * pow(Qvp,2)*Mt / (4*M_PI * pow( 2*ErGeV*Mt + pow(pList->mMed,2) ,2) ) 
	                     * ( 2*intConst - intInvEnuSq * ErGeV * Mt ) ); 
        }
        //axialvector
        case 4:
        {
            return convFactor * ( 
                       pList->SMinterference1 * GFERMI*Mt*Qa*Qap / (2*sqrt(2)*M_PI * ( 2*ErGeV*Mt + pow(pList->mMed,2) ) ) 
	                     * ( 2*intConst + intInvEnuSq * ErGeV * Mt ) 
	                 - pList->SMinterference1 * GFERMI*MN*Qv*Qap*ErGeV / (sqrt(2)*M_PI * ( 2*ErGeV*Mt + pow(pList->mMed,2) ) ) 
	                     * intInvEnu 
	                 + pList->SMinterference2 * pow(Qap,2)*Mt / (4*M_PI * pow( 2*ErGeV*Mt + pow(pList->mMed,2) ,2) ) 
	                     * ( 2*intConst + intInvEnuSq * ErGeV * Mt ) ); 
        }
        default:
        {
            return 1e-99;
        }
    }
    	
}

//full rate for a detector in events/kg/day/keV
double BSMrate(double ErKeV, paramList *pList, int detj)
{

    double rate = 0;
    double targetsPerKG;
    
    for(int i=0;i<pList->detectors[detj].nIso;i++)
	{   
        targetsPerKG = GeVperKG/(MN*pList->detectors[detj].isoA[i]); //how many targets per kg of detector
		
    	if(pList->nucScat)
    	{
    	    pList->Qa = pList->detectors[detj].isoSN[i]*(-0.427*-0.501163+0.842*0.506875) + pList->detectors[detj].isoSZ[i]*(-0.427*0.506875+0.842*-0.501163);	 
		    pList->Qv = ( (pList->detectors[detj].isoA[i] - pList->detectors[detj].isoZ[i]) - (1-4*0.2312*pList->detectors[detj].isoZ[i]) ) * ffactorSI( pList->detectors[detj].isoA[i], ErKeV);	 
		    pList->Qs = ((pList->detectors[detj].isoA[i] - pList->detectors[detj].isoZ[i]) * pList->qNs + pList->detectors[detj].isoZ[i] * pList->qPs ) * ffactorSI( pList->detectors[detj].isoA[i], ErKeV);
		    pList->Qvp = ((pList->detectors[detj].isoA[i]-pList->detectors[detj].isoZ[i]) * pList->qNv + pList->detectors[detj].isoZ[i] * pList->qPv) * ffactorSI( pList->detectors[detj].isoA[i], ErKeV);
		    pList->Qap = pList->detectors[detj].isoSN[i] * pList->qNa + pList->detectors[detj].isoSZ[i] * pList->qPa;

		    rate += targetsPerKG * pList->detectors[detj].isoFrac[i] * BSMrateN( ErKeV, pList, MN*pList->detectors[detj].isoA[i]);
	    }
	    if(pList->elecScat)
	    {
            pList->qV = 0.5+2*0.2312;
	        pList->qA = 0.5;
		    rate += pList->detectors[detj].isoZ[i] * targetsPerKG * pList->detectors[detj].isoFrac[i] * BSMrateE( ErKeV, pList, ME);
		}
	}	
	
	return rate;
}   

//template BSM rate
double aBSMrate(double ErKeV, paramList *pList, int detj)
{
	double rate = 0;
	paramList pListBSM = *pList;
   
	//nucleon axial charges
	pListBSM.qAp = 0;
	pListBSM.qAn = 0;
	
	//nucleon vector charges
	pListBSM.qVp = 2*pListBSM.qVu + pListBSM.qVd;
	pListBSM.qVn = pListBSM.qVu + 2*pListBSM.qVd;
	
	//weighted sum over different isotopes
	for(int i=0;i<pList->detectors[detj].nIso;i++)
	{
		//nuclei charges
		pListBSM.qA = pListBSM.qAn * pList->detectors[detj].isoSN[i] + pListBSM.qAp * pList->detectors[detj].isoSZ[i];	
		pListBSM.qV = ( pListBSM.qVn * (pList->detectors[detj].isoA[i] - pList->detectors[detj].isoZ[i]) + pListBSM.qVp * pList->detectors[detj].isoZ[i] ) * ffactorSI( pList->detectors[detj].isoA[i], ErKeV);	   
	   
		rate += pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListBSM, MN*pList->detectors[detj].isoA[i]);
	}
	   
	return rate; 
}

double E6rate(double ErKeV, paramList *pList, int detj)
{
	double rate = 0;
	paramList pListBSM = *pList;
   
    double Mzp = 1e3; //GeV)
    if(pList->BSM == 1)
    {
        pListBSM.qVu=0.191004;
        pListBSM.qVd=-0.351608-1939.0/pow(Mzp,2);
        pListBSM.qAu=-0.501163+969.5/pow(Mzp,2);
        pListBSM.qAd=0.506875-969.5/pow(Mzp,2);
    }
    else
    {
        pListBSM.qVu=0;
        pListBSM.qVd=0;
        pListBSM.qAu=0;
        pListBSM.qAd=0;
    }
	//nucleon axial charges
	pListBSM.qAp =  0.842*pListBSM.qAu-0.427*pListBSM.qAd;
	pListBSM.qAn = -0.427*pListBSM.qAu+0.842*pListBSM.qAd;
	
	//nucleon vector charges
	pListBSM.qVp = 2*pListBSM.qVu + pListBSM.qVd;
	pListBSM.qVn = pListBSM.qVu + 2*pListBSM.qVd;

	//weighted sum over different isotopes
	for(int i=0;i<pList->detectors[detj].nIso;i++)
	{
		//nuclei charges
		pListBSM.qA = pListBSM.qAn * pList->detectors[detj].isoSN[i] + pListBSM.qAp * pList->detectors[detj].isoSZ[i];	
		pListBSM.qV = ( pListBSM.qVn * (pList->detectors[detj].isoA[i] - pList->detectors[detj].isoZ[i]) + pListBSM.qVp * pList->detectors[detj].isoZ[i] ) * ffactorSI( pList->detectors[detj].isoA[i], ErKeV);	   
	   
		rate += pList->detectors[detj].isoFrac[i] * nuRate( ErKeV, &pListBSM, MN*pList->detectors[detj].isoA[i]);
	}

	return rate; 
	
}

double diffBSMrate(double ErkeV, paramList *pList, int detj, double signalNorm)						  
{   
    double rate=0;
    for(int i=0; i< pList->source.numFlux; i++)
    {
        if(pList->BSM==3 || pList->BSM==4)
        {
            rate += pList->source.nuFluxNorm[i] * ( signalNorm * gsl_spline_eval(pList->detectors[detj].signalBSM1[i], ErkeV, pList->detectors[detj].accelBSM1[i])
                + ((signalNorm > 0) - (signalNorm < 0) )*pow(signalNorm,2) * gsl_spline_eval(pList->detectors[detj].signalBSM2[i], ErkeV, pList->detectors[detj].accelBSM2[i]) );
        }
        else
        {
            rate += pList->source.nuFluxNorm[i] * signalNorm * gsl_spline_eval(pList->detectors[detj].signalBSM1[i], ErkeV, pList->detectors[detj].accelBSM1[i]);
        }
    }
    return rate;
}

double intBSMrate(double Er_min, double Er_max, paramList *pList, int detj, double signalNorm)						  
{   
    double rate=0;
    for(int i=0; i< pList->source.numFlux; i++)
    {
        if(pList->BSM==3 || pList->BSM==4)
        {
            rate +=  pList->source.nuFluxNorm[i] *( signalNorm * gsl_spline_eval_integ(pList->detectors[detj].signalBSM1[i], Er_min, Er_max, pList->detectors[detj].accelBSM1[i])
                + ((signalNorm > 0) - (signalNorm < 0) )*pow(signalNorm,2) * gsl_spline_eval_integ(pList->detectors[detj].signalBSM2[i], Er_min, Er_max, pList->detectors[detj].accelBSM2[i]) );
        }
        else
        {
            rate += pList->source.nuFluxNorm[i] * signalNorm * gsl_spline_eval_integ(pList->detectors[detj].signalBSM1[i], Er_min, Er_max, pList->detectors[detj].accelBSM1[i]);
        }
    }
    return rate;
}
