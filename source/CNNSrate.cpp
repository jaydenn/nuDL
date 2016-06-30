#include <iostream>
#include <cmath>
#include "DEIntegrator.h"
#include "parameterStruct.h"
#include "formFactorSI.h"
#include "nuFlux.h"
#ifndef GSL_INTERP_H
    #include <gsl/gsl_interp.h>
#endif
#ifndef GSL_INTEGRATION_H
    #include <gsl/gsl_integration.h>
#endif 
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif  

double FIRSTEVAL = 0;
gsl_integration_workspace * W;
const double GFERMI = 1.1664e-5; //GeV^-2
const double HBAR   = 6.5821e-25; //GeV*s
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
    return nuFlux(EnuGeV, pList) /  EnuGeV ;
}

double EnuIntegrand2(double EnuGeV, void *pars)
{
    paramList *pList = (paramList*)pars;
    return nuFlux(EnuGeV, pList) / ( EnuGeV*EnuGeV ) ;
}

double CNSintegral(double ErGeV,  paramList *pList, int detj, int EnuPow)
{
    int key = 2;
    int limit = 1000;
    double integral,absErr;
    
    pList->A = pList->detectors[detj].AM;
    pList->Er = ErGeV;
            
    if(FIRSTEVAL==0)
    {
        W = gsl_integration_workspace_alloc (1000);
        FIRSTEVAL=1;
    }
    
    if(EnuPow==0)
        pList->F.function = &EnuIntegrand0;
    else if(EnuPow==-1)
        pList->F.function = &EnuIntegrand1;
    else if(EnuPow==-2)
        pList->F.function = &EnuIntegrand2;
             
    pList->F.params = pList; //yeah, that's not weird..
    
    double EnuMinGeV = sqrt( ErGeV * (MN*pList->detectors[detj].AM) / 2);

    gsl_integration_qag(&(pList->F), EnuMinGeV, pList->EnuMax, 1e-22, 1e-4, limit, 3, W, &integral, &absErr); //check normalization

    return integral;    
}

//returns CNNS rate per kg/year/keV
double CNNSrate(double ErKeV, paramList *pList, int detj) 						  
{

    double ErGeV = ErKeV/GeVtoKeV;
    double secsPerYear = 365.25*24*60*60; 
    double atomsPerKG =  GeVperKG / ( pList->detectors[detj].AM * MN); 

    double intConst    = CNSintegral( ErGeV, pList, detj, 0);
    double intInvEnu   = CNSintegral( ErGeV, pList, detj, -1);
    double intInvEnuSq = CNSintegral( ErGeV, pList, detj, -2);

    return pow(GFERMI,2) / ( 2 * M_PI ) * (MN*pList->detectors[detj].AM) / GeVtoKeV * atomsPerKG * secsPerYear
            * ( 
                  intConst    * 2*( pow(pList->qA,2) + pow(pList->qV,2) )
                - intInvEnu   * 2*ErGeV*pow(pList->qA-pList->qV,2) 
                + intInvEnuSq * ( ErGeV*ErGeV*pow(pList->qA - pList->qV,2) + ErGeV*MN*pList->detectors[detj].AM*(pow(pList->qA,2) - pow(pList->qV,2)) )
              ); 
    
}

double BSMrate(double ErKeV, paramList *pList, int detj)
{
    double rate = 0;
    paramList pListBSM = *pList;
   
    //nucleon axial charges
    pListBSM.qAp = 0;
    pListBSM.qAn = 0;
    
    //nucleon vector charges
    pListBSM.qVp = 0;
    pListBSM.qVn = 0;
    
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
        pListSM.qA = 0.5 * (pList->detectors[detj].isoSN[i] - pList->detectors[detj].isoSZ[i]);  
        pListSM.qV = 0.5 * (pList->detectors[detj].isoA[i] - 1.04*pList->detectors[detj].isoZ[i])* ffactorSI( pList->detectors[detj].isoA[i], ErKeV);    
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
        return CNNSrate( Er, pList, detj);
    }
};

double intCNNSrate(double Er_min, double Er_max, paramList *pList, int detj)
{
    ErIntegral erInt; 		       //create instance for recoil energy integration
    erInt.pList = pList;
    erInt.detj = detj;

    return DEIntegrator<ErIntegral>::Integrate(erInt,Er_min,Er_max,1e-4);
}
