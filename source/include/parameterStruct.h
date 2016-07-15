//definition of parameter structs
#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>
#include <string.h>
#ifndef IOSTREAM
	#include <iostream>
#endif
#ifndef DETECTORSTRUCT_H
	#include "detectorStruct.h"
#endif	
#define PARAMETERSTRUCT_H

struct paramList {
	
	//root directory for file output
	char root[50];
	
	//arrays for reactor flux data
	gsl_interp *nuFluxInterp;
	gsl_interp_accel *nuFluxAccel;
	double nuFlux, nuFluxUn;
	double nuFluxNorm, signalNorm;
	double Er, EnuMax;
    
	//setup for rate integration
	gsl_function F;
	double (*rateFunc)(double, double, paramList *, int);
     
	double A; //is this used?
	double qA, qAn, qAp, qAu, qAd;
	double qV, qVn, qVp, qVu, qVd;
	
	//fiducial coupling
	double C;
	
	//nucleon couplings
	double Qs;  //scalar nuclear
	double qPs; //scalar
	double qNs; //scalar
	double Qv;  //vector nuclear
    double qPv; //vector
    double qNv; //vector
    double Qa;  //axial nuclear
    double qPa; //axial
    double qNa; //axial
    
    //neutrino couplings
    double gNuV;  //vector neutrino
    double gNuS;  //scalar neutrino
        
    //electron couplings	
	double gEs;   //scalar electron
    double gEp;   //pseudoscalar electron
    double gEv;  //vector electron
    double gEa;  //axial-vector electron
    
	int BSM;
	int nucScat;
	int elecScat;
    int mediator;
    double mMed;
    
    int asimov;
	int ndet, detj;
	detector detectors[10];
	
	void printPars()
	{
   		std::cout << "Conudl configuration:" << std::endl;		
    	std::cout << "  root: " << root << std::endl;
    	std::cout << "  BSM: " << BSM << std::endl;
    	std::cout << "  asimov: " << asimov << std::endl;
    	std::cout << "  flux: " << nuFlux << " +/- " << nuFluxUn*100 << "%" << std::endl;
		std::cout << ndet << " detectors:" << std::endl;
		for(int i=0;i<ndet;i++)
			detectors[i].printDetSpecs();
	}
	
	paramList()
	{
	    Qs=qNs=qPs=qNv=qPv=qNa=qPa=gNuS=gNuV=gEs=gEp=gEv=gEa=C=0;
		ndet=0;
		nuFluxNorm=1;
		nuFluxUn=1e-99;
        signalNorm=1;
		EnuMax=0;
        asimov=1;
        elecScat=0;
        nucScat=0;
        mediator=0;
	}
	
};
	
