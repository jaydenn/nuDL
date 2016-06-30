#include <iostream>
#include <string.h>
#include <cmath>
#include "DEIntegrator.h"
#include "detectorFunctions.h"
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

//returns integrated total # background events per tonne/year for bg type, with recoil  Er_min < Er < Er_max
double intBgRate(detector det, double Er_min, double Er_max)  						  
{   
    return 0; //gsl_spline_eval_integ(det.background, Er_min, Er_max, det.accel);
}

double diffBgRate(detector det, double Er)
{
    return 0; //gsl_spline_eval(det.background, Er, det.accel);
}

//only need to calculate background once, store in a table for interpolation, stored as N/t/year/keV
int InitializeBackground(detector *det)
{
    //get values of bg at relevant energies
    double Er[500];
    double Bg[500];
    for(int i=0; i<500; i++)
    {
        Er[i] = det->ErL + (double)i*(det->ErU-det->ErL)/499;
        Bg[i] = 0; //detBackground(Er[i], det->bg) * 1000.0 * 365.24;
    }
    
    //create gsl interpolation object
    return gsl_spline_init(det->background,Er,Bg,500);

}

int newDetector(detector *det, char *name, double exp, int ndet)
{
    det->exposure = exp;

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
    
    sprintf( det->name, "%s", &(name[1]));
    int nIso = det->nIso = 0;
    
    while(temp[0]!='-')
    {
        err = fscanf(detsINI,"%s",temp);
                     
        if(strcmp(temp,"AM")==0)  
            err=fscanf(detsINI,"%lf",&(det->AM)); 
        if(strcmp(temp,"Er")==0)  
            err=fscanf(detsINI,"%lf-%lf",&(det->ErL),&(det->ErU)); 
        if(strcmp(temp,"bg")==0)  
            err=fscanf(detsINI,"%d",&(det->bg));
        if(strcmp(temp,"eff")==0) 
            err=fscanf(detsINI,"%d",&(det->eff));         
        if(strcmp(temp,"res")==0) 
            err=fscanf(detsINI,"%d",&(det->res));
    }

    ret = fgets(temp,200,detsINI);
    ret = fgets(temp,200,detsINI);
    ret = fgets(temp,200,detsINI);
    ret = fgets(temp,200,detsINI);
    while(!feof(detsINI) && temp[0]!='_')
    {       
        sscanf(temp,"%d %d %lf %lf %lf",&(det->isoZ[nIso]),&(det->isoA[nIso]),&(det->isoFrac[nIso]),&(det->isoSZ[nIso]),&(det->isoSN[nIso])); 
        nIso++;
        ret = fgets(temp,200,detsINI);         
    }
    det->nIso = nIso;
             
    fclose(detsINI);
    
    if(InitializeBackground(det))
    { 
        std::cout << "problem initializing background of detector " << name << std::endl;
        return 1;
    }

    return 0;
    
}


