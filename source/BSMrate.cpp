#include <iostream>
#include <cmath>
#ifndef GSL_INTERP_H
	#include <gsl/gsl_interp.h>
#endif
#include "formFactorSI.h"
#include "nuRate.h"
#include "nuFlux.h"
#include "sterileRate.h"
#include "detectorFunctions.h"
#ifndef DETECTORSTRUCT_H
	#include "detectorStruct.h"
#endif
#ifndef PARAMETERSTRUCT_H
	#include "parameterStruct.h"
#endif	

#include "physicalConstants.h"

//returns simplified model light mediator rate per electron/day/keV
double BSMrateE(double ErKeV, paramList *pList, double Mt, int fluxj)						  
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
        
    intConst	= fluxIntegral( ErGeV, pList, ME,  0, fluxj);  //units GeV^2/s
    intInvEnuSq = fluxIntegral( ErGeV, pList, ME, -2, fluxj); //units 1/s
    	            
    switch( pList->BSM )
    {
        //scalar
        case 1:
        {
            return convFactor * pow(gEs*gNuS,2) * ErGeV * pow(ME,2) / ( 4*M_PI )              
                * intInvEnuSq / pow( 2*ErGeV*ME + pow(pList->mMed,2) ,2);
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
                + pList->SMinterference2 * ( pow(gEv*gNuV,2) * ME / ( 2 * M_PI ) 
                * intConst / pow( 2*ErGeV*ME + pow(pList->mMed,2) ,2) 
                + 0*pow(gNuV*gNuV,2) * ME  * pow(ErGeV,2) / ( 4 * M_PI ) 
                * intInvEnuSq / pow( 2*ErGeV*ME + pow(pList->mMed,2) ,2)) );   
        }
        //axialvector
        case 4:
        {
            return convFactor * ( pList->SMinterference1 * sqrt(2)*GFERMI*ME*ga*gNuV*gEa /  M_PI 
                * intConst / ( 2*ErGeV*ME + pow(pList->mMed,2) )
                + pList->SMinterference2 * ( pow(gEa*gNuV,2)* ME / ( 2 * M_PI ) 
                * intConst / pow( 2*ErGeV*ME + pow(pList->mMed,2) ,2)
                + 0*pow(gNuV*gNuV,2)* ME * pow(ErGeV,2) / ( 4 * M_PI ) 
                * intInvEnuSq / pow( 2*ErGeV*ME + pow(pList->mMed,2) ,2) ) );   
        }
        default:
        {
            return 1e-99;
        }
    }
}

//returns simplified model light mediator rate per nuclei/day/keV
double BSMrateN(double ErKeV, paramList *pList, double Mt, int fluxj)						  
{

	double ErGeV = ErKeV/GeVtoKeV;
	//double secsPerDay = 24*60*60; 
	
    double convFactor =  1/GeVtoKeV * secsPerDay; //units GeV*s/keV/kg/day
    
    //nuclear couplings
    double Qv = pList->Qv; //SM vector
    double Qa = pList->Qa; //SM axial coupling
    
    //makes calculations more readable
    double Qs = pList->Qs;  //scalar
    double Qvp = pList->Qvp; //vector
    double Qap = pList->Qap; //axial 

    double intConst, intInvEnu, intInvEnuSq;
    
    intConst	= fluxIntegral( ErGeV, pList, Mt,  0, fluxj); //units GeV^2/s
    intInvEnu   = fluxIntegral( ErGeV, pList, Mt, -1, fluxj); //units GeV/s
    intInvEnuSq = fluxIntegral( ErGeV, pList, Mt, -2, fluxj); //units 1/s

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
	                 - pList->SMinterference1 * GFERMI*Mt*Qv*Qap*ErGeV / (sqrt(2)*M_PI * ( 2*ErGeV*Mt + pow(pList->mMed,2) ) ) 
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
double BSMrate(double ErKeV, paramList *pList, int detj, int fluxj)
{

    double rate = 0;
    double targetsPerKG;
    if(pList->BSM==0)
        return 0;
    if(pList->BSM == 6)
        return sterileRate(ErKeV, pList, detj, fluxj);
    
    for(int i=0;i<pList->detectors[detj].nIso;i++)
	{   
        targetsPerKG = pList->detectors[detj].isoFrac[i] *GeVperKG/(AMU*pList->detectors[detj].isoA[i]); //how many targets per kg of detector
    	if(pList->nucScat)
    	{
    	    pList->Qa = 4.0/3.0 * (pList->detectors[detj].isoJN[i]+1) / pList->detectors[detj].isoJN[i] * ( pList->detectors[detj].isoSN[i]*GAN + pList->detectors[detj].isoSZ[i]*GAP );	 
		    pList->Qv = ( GVN * (pList->detectors[detj].isoA[i] - pList->detectors[detj].isoZ[i]) + GVP * pList->detectors[detj].isoZ[i] )* ffactorSI( pList->detectors[detj].isoA[i], ErKeV);
		    pList->Qvp = ((pList->detectors[detj].isoA[i]-pList->detectors[detj].isoZ[i]) * pList->qNv + pList->detectors[detj].isoZ[i] * pList->qPv) * ffactorSI( pList->detectors[detj].isoA[i], ErKeV);
		    pList->Qap = 4.0/3.0 * (pList->detectors[detj].isoJN[i]+1) / pList->detectors[detj].isoJN[i] * ( pList->detectors[detj].isoSN[i] * pList->qNa + pList->detectors[detj].isoSZ[i] * pList->qPa );
            
		    rate += targetsPerKG * BSMrateN( ErKeV, pList, MN*pList->detectors[detj].isoA[i], fluxj);
	    }
	    if(pList->elecScat)
	    {
	        int Ne=0;
	        while(pList->detectors[detj].ionization[i][Ne] > ErKeV && Ne < pList->detectors[detj].isoZ[i]) 
	            Ne++;
            if(pList->source.isSolar[fluxj] == 1)
	        {
		        pList->qA = 0.5;
		        pList->qV = 2*SSW+0.5;
		        rate += pList->source.survProb[fluxj] * ((double) pList->detectors[detj].isoZ[i] - Ne) * targetsPerKG  * BSMrateE( ErKeV, pList, ME, fluxj);
		        
		        pList->qA = -0.5;
		        pList->qV = 2*SSW-0.5;
                rate += (1-pList->source.survProb[fluxj]) * ((double) pList->detectors[detj].isoZ[i] - Ne) * targetsPerKG * BSMrateE( ErKeV, pList, ME, fluxj);		     
	        }
		    else
		    {
		        //only include if doing reactor neutrinos (anti-neutrino channel)
		        if(pList->source.isSolar[0] == 0)
		        {
		            pList->qA = -0.5;
		            pList->qV = 2*SSW+0.5;
                    rate += ((double) pList->detectors[detj].isoZ[i] - Ne) * targetsPerKG * BSMrateE( ErKeV, pList, ME, fluxj);
                }
		    }
	
		}
	}	

	return rate*detEff(ErKeV,pList->detectors[detj].eff);
}


double diffBSMrate(double ErkeV, paramList *pList, int detj, double signalNorm)						  
{   
    double rate=0;
    for(int i=0; i< pList->source.numFlux; i++)
    {
        if(pList->BSM==3 || pList->BSM==4)
            rate += pList->source.nuFluxNorm[i] * ( signalNorm * gsl_spline_eval(pList->detectors[detj].signalBSM1[i], ErkeV, pList->detectors[detj].accelBSM1[i])
                + ((signalNorm > 0) - (signalNorm < 0) )*pow(signalNorm,2) * gsl_spline_eval(pList->detectors[detj].signalBSM2[i], ErkeV, pList->detectors[detj].accelBSM2[i]) );
        else
            rate += pList->source.nuFluxNorm[i] * signalNorm * gsl_spline_eval(pList->detectors[detj].signalBSM1[i], ErkeV, pList->detectors[detj].accelBSM1[i]);
    }
    return rate;
}

double intBSMrate(double Er_min, double Er_max, paramList *pList, int detj, double signalNorm)						  
{   
    double rate=0;
    for(int i=0; i< pList->source.numFlux; i++)
    {
        if(pList->BSM==3 || pList->BSM==4)
            rate +=  pList->source.nuFluxNorm[i] * ( signalNorm * gsl_spline_eval_integ(pList->detectors[detj].signalBSM1[i], Er_min, Er_max, pList->detectors[detj].accelBSM1[i])
                + ((signalNorm > 0) - (signalNorm < 0) )*pow(signalNorm,2) * gsl_spline_eval_integ(pList->detectors[detj].signalBSM2[i], Er_min, Er_max, pList->detectors[detj].accelBSM2[i]) );
        else
            rate += pList->source.nuFluxNorm[i] * signalNorm * gsl_spline_eval_integ(pList->detectors[detj].signalBSM1[i], Er_min, Er_max, pList->detectors[detj].accelBSM1[i]);
    }
    
    return rate;
}
