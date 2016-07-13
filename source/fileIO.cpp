#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string.h> 
#include <stdlib.h>
#include "detectorFunctions.h"
#include "nuFlux.h"
#include "SMrate.h"
#include "BSMrate.h"
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif
#ifndef DETECTORSTRUCT_H
	#include "detectorStruct.h"
#endif	

//gets sampling parameters from file
int readConfigFile(paramList *pL, char *filename) 
{

    FILE* input;
    input = fopen(filename,"r");
    if(input==NULL) 
    {
        printf("unable to open parameter file: %s\n",filename);
        return -1;	
    }
  
    int mode;
    char *ret;
    char temp[400];
    char root[50];
    ret = fgets(temp,200,input);

    //Mode switch
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&mode);
    
    // root for output files
    ret = fgets(temp,200,input);
    sscanf(temp,"%s %*s",root);
    sprintf(pL->root, "%s", root);		 
  
    //include nuclear scattering?
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&(pL->nucScat));
    
    //include electron scattering?
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&(pL->elecScat));
    
    if( (pL->nucScat == 0 && pL->elecScat == 0) || (pL->nucScat == 1 && pL->elecScat == 1))
    {
        std::cout << "Must choose one of electron or nuclear scattering" << std::endl;
        return -1;
    }
    
    //which BSM model to consider
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&(pL->BSM));
    if(pL->BSM==0)
    {
        std::cout << "Must choose a BSM type" << std::endl;
        return -1;
    }
    else
        pL->rateFunc = *intBSMrate;
    
    //mediator mass
    ret = fgets(temp,200,input);
    sscanf(temp,"%lf",&(pL->mMed));
    
    //reactor flux
    ret = fgets(temp,200,input);
    ret = fgets(temp,200,input);
    sscanf(temp,"%lf %*s %lf",&(pL->nuFlux),&(pL->nuFluxUn));
    pL->nuFluxUn /= pL->nuFlux; //want fractional uncertainty
    
    //initialize reactor flux
	double distance = 1.0;
    pL->nuFlux /= pow(distance,2); 
    int err = initFlux(pL);
    if (err == 1)
    {
        std::cout << "Reactor flux data is not properly normalized" << std::endl;
        return -1;
    }
    
    //Detector setup
    char name[20];
    double exp;
    ret = fgets(temp,200,input);
    ret = fgets(temp,200,input);

    while(temp[0]=='#')
    {
        sscanf(temp,"%s %lf", name, &exp);
        
	    if(newDetector(pL, name, exp)) { std::cout << "Could not create detector, exiting" << std::endl; return -1; }
        std::cout << "Using detector " << name << std::endl;
        ret = fgets(temp,200,input);
    }    
        
    //asimov or random sim?
    ret = fgets(temp,200,input);
    ret = fgets(temp,200,input);
    sscanf(temp,"%d",&(pL->asimov));
    
    fclose(input); 

    return mode;
}

