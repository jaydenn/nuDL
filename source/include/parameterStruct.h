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
#ifndef DETECTORFUNCTIONS_H
	#include "detectorFunctions.h"
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
	double qA, qAn, qAp;
	double qV, qVn, qVp;
	
	int BSM;
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
		std::cout << "Detectors:" << std::endl;
		for(int i=0;i<ndet;i++)
			detectors[i].printDetSpecs();
	}
	
	paramList()
	{
		ndet=0;
		nuFluxNorm=1;
		nuFluxUn=1e-99;
        signalNorm=1;
		EnuMax=0;
        asimov=1;
	}
	
	int newDetector(char *name, double exp)
	{
		if(ndet==10)
		{
			std::cout << "already at max number of detectors (10)" << std::endl;
			return 1;
		}
			
		detector newDet;
		
		newDet.exposure = exp;
		
		//read in detector configuration
		    FILE *detsINI;
		    detsINI = fopen("detectors.ini","r");
		    if(detsINI==NULL)
		    {
			    printf("unable to open detectors.ini\n");
			    return 1;
		    }

		    char temp[200];
		    char *ret;
		    int err;
		
		    ret = fgets(temp,200,detsINI);
		
		    while(strcmp(temp,name)!=0)
		    {
			    err = fscanf(detsINI,"%s",temp);
			
			    if(feof(detsINI))
			    {
				    printf("detector '%s' not found\n",name); 
				    fclose(detsINI);
				    return 1;
			    }
		    }
		
		    sprintf( newDet.name, "%s", &(name[1]));
		
		    while(temp[0]!='-')
		    {
			    err = fscanf(detsINI,"%s",temp);
						     
			    if(strcmp(temp,"AM")==0)  
				    err=fscanf(detsINI,"%lf",&(newDet.AM)); 
			    if(strcmp(temp,"Er")==0)  
				    err=fscanf(detsINI,"%lf-%lf",&(newDet.ErL),&(newDet.ErU)); 
			    if(strcmp(temp,"bg")==0)  
				    err=fscanf(detsINI,"%d",&(newDet.bg));
				if(strcmp(temp,"bgUn")==0)  
				    err=fscanf(detsINI,"%lf",&(newDet.BgUn));
			    if(strcmp(temp,"eff")==0) 
				    err=fscanf(detsINI,"%d",&(newDet.eff));			
			    if(strcmp(temp,"res")==0) 
				    err=fscanf(detsINI,"%d",&(newDet.res));
		    }

		    ret = fgets(temp,200,detsINI);
		    ret = fgets(temp,200,detsINI);
		    ret = fgets(temp,200,detsINI);
		    ret = fgets(temp,200,detsINI);
		
		    while(!feof(detsINI) && temp[0]!='_')
		    {	
			    if(newDet.nIso==10)
			    {
				    std::cout << "already at max number of isotopes (10)" << std::endl;
				    break;
			    }
			    sscanf(temp,"%d %d %lf %lf %lf",&(newDet.isoZ[newDet.nIso]),&(newDet.isoA[newDet.nIso]),&(newDet.isoFrac[newDet.nIso]),&(newDet.isoSZ[newDet.nIso]),&(newDet.isoSN[newDet.nIso])); 
			
			    newDet.nIso++;
			    ret = fgets(temp,200,detsINI);		   
		    }
		    //finished reading in det data		 
		fclose(detsINI);
		
		//only need to calculate background and SM signal once, store in a table for interpolation, stored as events/kg/day/keV
		//get values of bg at relevant energies
		double Er[1000];
		double Bg[1000];
		for(int i=0; i<1000; i++)
		{
			Er[i] = newDet.ErL + (double)i*(newDet.ErU-newDet.ErL)/900;
			Bg[i] = detBackground(Er[i], newDet.bg);
		}
	
		//create gsl interpolation object
		gsl_spline_init(newDet.background,Er,Bg,1000);
		
		//add new detector to the array
		detectors[ndet] = newDet;
		ndet++;
		
		return 0;
		
	}
};
	
