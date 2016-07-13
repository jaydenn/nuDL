#include <iostream>
#include <cstdio>
#include <cmath>
#ifndef DETECTORSTRUCT_H
	#include "detectorStruct.h"
#endif	
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif
#include "nuRate.h"
#include "SMrate.h"
#include "BSMrate.h"

double detEff(double Er, int type)
{
    switch( type ) 
    {
        case 0: 
            return 1;
        case 1: 
            return .5;
        default:
            printf("invalid detector efficiency\n"); 
            return NAN;
    }
}

//return background in events/kg/day/keV
double detBackground(double Er, paramList *pList, int detj)
{
    
    switch( pList->detectors[detj].bg ) 
    {
        case 0: 
            return 1e-99;
        case 1: 
            return 100.0;
	    case 2: 
            return 10.0;
        default:
            printf("invalid detector background\n"); 
            return NAN; 
    }
}

double detRes(double Er, int type)
{
    switch( type ) 
    {
        case 0: 
            return 0;
        case 1: 
            return 0;
        case 2: 
            return 0;
        case 3: 
            return 0;
        case 4: 
            return 0;
        default:
            printf("invalid detector resolution\n"); 
            return NAN; 
    }
}

int newDetector(paramList *pList, char *name, double exp)
{
		if(pList->ndet==10)
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
		
		//add new detector to the array
		pList->detectors[pList->ndet] = newDet;
		
		//only need to calculate background and SM signal once, store in a table for interpolation, stored as events/kg/day/keV
		//get values of bg at relevant energies
		std::cout << "Initializing rates:" << std::endl << "SM rate..." << std::endl; 
		rateInit( pList, pList->ndet,        &SMrate,   pList->detectors[pList->ndet].signalSM);
		if(pList->BSM!=0)
		{   
		    std::cout << "BSM rate..." << std::endl; 
		    rateInit( pList, pList->ndet,       &BSMrate,  pList->detectors[pList->ndet].signalBSM);
        }
        std::cout << "Background rate..." << std::endl; 
		rateInit( pList, pList->ndet, &detBackground, pList->detectors[pList->ndet].background);
		std::cout << "done." << std::endl; 
		
		pList->ndet++;
		
		return 0;
		
}

//returns integrated total # background events per tonne/year for bg type, with recoil  Er_min < Er < Er_max
double intBgRate(detector det, double Er_min, double Er_max)  						  
{   
    return gsl_spline_eval_integ(det.background, Er_min, Er_max, det.accelBg);
}

double diffBgRate(detector det, double Er)
{
    return gsl_spline_eval(det.background, Er, det.accelBg);
}
