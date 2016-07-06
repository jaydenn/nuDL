#include <iostream>
#include <cmath>
#ifndef GSL_INTERP_H
	#include <gsl/gsl_interp.h>
#endif
#ifndef GSL_INTEGRATION_H
	#include <gsl/gsl_integration.h>
#endif 
#include "DEIntegrator.h"
#include "formFactorSI.h"
#include "nuFlux.h"
#ifndef DETECTORSTRUCT_H
	#include "detectorStruct.h"
#endif
#ifndef PARAMETERSTRUCT_H
	#include "parameterStruct.h"
#endif	

double FIRSTEVALCNS = 0;
double FIRSTEVALSM = 0;
gsl_integration_workspace * W;
const double GFERMI = 1.1664e-5; //GeV^-2
const double HBAR	= 6.5821e-25; //GeV*s
const double GeVtoKeV = 1e6;
const double GeVperKG = 5.6094e26;
const double MN = 0.9383; //mass of nucleon in GeV

double EnuIntegrand0(double EnuGeV, void *pars)
{
	paramList *pList = (paramList*)pars;
	return nuFlux(EnuGeV, pList);
}

double EnuIntegrand1(double EnuGeV, void *pars)
{
	paramList *pList = (paramList*)pars;
	return nuFlux(EnuGeV, pList) /	EnuGeV ;
}

double EnuIntegrand2(double EnuGeV, void *pars)
{
	paramList *pList = (paramList*)pars;
	return nuFlux(EnuGeV, pList) / ( EnuGeV*EnuGeV ) ;
}

double CNSintegral(double ErGeV,  paramList *pList, int detj, int EnuPow)
{
	int key = 2;
	int limit = 2000;
	double integral,absErr,tol;
	
	pList->A = pList->detectors[detj].AM;
	pList->Er = ErGeV;
			
	if(FIRSTEVALCNS==0)
	{
		W = gsl_integration_workspace_alloc (2000);
		FIRSTEVALCNS=1;
	}
	
	if(EnuPow==0)
	{
		pList->F.function = &EnuIntegrand0;
		tol=1e-22;
	}
	else if(EnuPow==-1)
	{
		pList->F.function = &EnuIntegrand1;
		tol=1e-20;
	}
	else if(EnuPow==-2)
	{
		pList->F.function = &EnuIntegrand2;
		tol=1e-18;
	}
			 
	pList->F.params = pList; //yeah, that's not weird..
	
	double EnuMinGeV = sqrt( ErGeV * (MN*pList->detectors[detj].AM) / 2);

	gsl_integration_qag(&(pList->F), EnuMinGeV, pList->EnuMax, tol, 1e-3, limit, 2, W, &integral, &absErr); 
//std::cout << ErGeV  << " " << EnuIntegrand2(ErGeV,(void *)pList) << " " << EnuIntegrand1(ErGeV,(void *)pList) << " " << EnuIntegrand0(ErGeV,(void *)pList) << std::endl;
	return integral;	
}

//returns CNNS rate per kg/day/keV
double CNNSrate(double ErKeV, paramList *pList, int detj)						  
{

	double ErGeV = ErKeV/GeVtoKeV;
	double secsPerDay = 24*60*60; 
	double atomsPerKG =	 GeVperKG / ( pList->detectors[detj].AM * MN); 

	double intConst	   = CNSintegral( ErGeV, pList, detj, 0);
	double intInvEnu   = CNSintegral( ErGeV, pList, detj, -1);
	double intInvEnuSq = CNSintegral( ErGeV, pList, detj, -2);

	return pow(GFERMI,2) / ( 2 * M_PI ) * (MN*pList->detectors[detj].AM) / GeVtoKeV * atomsPerKG * secsPerDay
			* ( 
				  intConst	  * 2*( pow(pList->qA,2) + pow(pList->qV,2) )
				- intInvEnu	  * 2*ErGeV*pow(pList->qA-pList->qV,2) 
				+ intInvEnuSq * ( ErGeV*ErGeV*pow(pList->qA - pList->qV,2) + ErGeV*MN*pList->detectors[detj].AM*(pow(pList->qA,2) - pow(pList->qV,2)) )
			  ); 
	
}

double BSMrate(double ErKeV, paramList *pList, int detj)
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
	   
		rate += pList->detectors[detj].isoFrac[i] * CNNSrate( ErKeV, &pListBSM, detj);
	}

	return rate; 
	
}

//returns SM rate per kg/year/keV
double SMrate(double ErKeV, paramList *pList, int detj)							  
{
    
	double rate = 0;
	paramList pListSM = *pList;
	
	for(int i=0;i<pList->detectors[detj].nIso;i++)
	{
		pListSM.qA = pList->detectors[detj].isoSN[i]*(-0.427*-0.501163+0.842*0.506875) + pList->detectors[detj].isoSZ[i]*(-0.427*0.506875+0.842*-0.501163);	 
		pListSM.qV = (- 0.512213 * (pList->detectors[detj].isoA[i] - pList->detectors[detj].isoZ[i]) + .0304*pList->detectors[detj].isoZ[i] )* ffactorSI( pList->detectors[detj].isoA[i], ErKeV);	 
		rate += pList->detectors[detj].isoFrac[i] * CNNSrate( ErKeV, &pListSM, detj);
	}

	return rate; 
}

class ErIntegral
{
public:
	paramList *pList;
	int detj;
	double operator()(double Er) const
	{
		return SMrate( Er, pList, detj);
	}
};

double intCNNSrate(double Er_min, double Er_max, paramList *pList, int detj)
{
	ErIntegral erInt;			   //create instance for recoil energy integration
	erInt.pList = pList;
	erInt.detj = detj;

	return DEIntegrator<ErIntegral>::Integrate(erInt,Er_min,Er_max,1e-4);
}

void rateInit( paramList *pList, int detj, double (*rateFunc)(double, paramList *, int), gsl_spline *rateSpline)	
{
    double ErkeV[1000];
    double rate[1000];
    for(int i=0; i<1000; i++)
    {
        ErkeV[i] = pList->detectors[detj].ErL + (double)i*(pList->detectors[detj].ErU-pList->detectors[detj].ErL)/900;
        rate[i] = rateFunc( (double)ErkeV[i], pList, detj);	
    }

    //create gsl interpolation object
    gsl_spline_init(rateSpline,ErkeV,rate,1000);
}

double diffSMrate(double ErkeV, paramList *pList, int detj)						  
{   

    if( ErkeV < (pList->detectors[detj].ErL + (double)999*(pList->detectors[detj].ErU-pList->detectors[detj].ErL)/900) )
        return gsl_spline_eval(pList->detectors[detj].signalSM, ErkeV, pList->detectors[detj].accelSM);
    else
        return 1e-99;
}

double intSMrate(double Er_min, double Er_max, paramList *pList, int detj)						  
{   
    return gsl_spline_eval_integ(pList->detectors[detj].signalSM, Er_min, Er_max, pList->detectors[detj].accelSM);
}

double diffBSMrate(double ErkeV, paramList *pList, int detj)						  
{   

    if( ErkeV < (pList->detectors[detj].ErL + (double)999*(pList->detectors[detj].ErU-pList->detectors[detj].ErL)/900) )
        return gsl_spline_eval(pList->detectors[detj].signalBSM, ErkeV, pList->detectors[detj].accelBSM);
    else
        return 1e-99;
}

double intBSMrate(double Er_min, double Er_max, paramList *pList, int detj)						  
{   
    return gsl_spline_eval_integ(pList->detectors[detj].signalBSM, Er_min, Er_max, pList->detectors[detj].accelBSM);
}
