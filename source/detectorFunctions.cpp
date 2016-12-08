#include <iostream>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <string>
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
double detBackground(double Er, paramList *pList, int detj, int fluxj)
{
    
    switch( pList->detectors[detj].bg ) 
    {
        case 0: 
            return 1e-99;
        case 1: 
            return 100.0;
	    case 2: 
            return 10.0;
        case 3: 
            return 1.0;
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

		pList->detectors[pList->ndet].exposure = exp;

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
	
	    sprintf( pList->detectors[pList->ndet].name, "%s", &(name[1]));
	
	    while(temp[0]!='-')
	    {
		    err = fscanf(detsINI,"%s",temp);
 
		    if(strcmp(temp,"Er")==0)  
			    err=fscanf(detsINI,"%lf-%lf",&(pList->detectors[pList->ndet].ErL),&(pList->detectors[pList->ndet].ErU)); 
		    if(strcmp(temp,"bg")==0)  
			    err=fscanf(detsINI,"%d",&(pList->detectors[pList->ndet].bg));
			if(strcmp(temp,"bgUn")==0)  
			    err=fscanf(detsINI,"%lf",&(pList->detectors[pList->ndet].BgUn));
		    if(strcmp(temp,"eff")==0) 
			    err=fscanf(detsINI,"%d",&(pList->detectors[pList->ndet].eff));			
		    if(strcmp(temp,"res")==0) 
			    err=fscanf(detsINI,"%d",&(pList->detectors[pList->ndet].res));
	    }
        
        //optimize ROI
        if (pList->nucScat == 1 && pList->detectors[pList->ndet].ErU > 8)
        {
            std::cout << "decreasing ROI to increase SNR\n";
            pList->detectors[pList->ndet].ErU = 12;
        }
        
	    ret = fgets(temp,200,detsINI);
	    ret = fgets(temp,200,detsINI);
	    ret = fgets(temp,200,detsINI);
	    ret = fgets(temp,200,detsINI);
	
	    while(!feof(detsINI) && temp[0]!='-')
	    {	
		    if(pList->detectors[pList->ndet].nIso==10)
		    {
			    std::cout << "already at max number of isotopes (10)" << std::endl;
			    break;
		    }
		    sscanf(temp,"%d %d %lf %lf %lf",&(pList->detectors[pList->ndet].isoZ[pList->detectors[pList->ndet].nIso]),&(pList->detectors[pList->ndet].isoA[pList->detectors[pList->ndet].nIso]),&(pList->detectors[pList->ndet].isoFrac[pList->detectors[pList->ndet].nIso]),&(pList->detectors[pList->ndet].isoSZ[pList->detectors[pList->ndet].nIso]),&(pList->detectors[pList->ndet].isoSN[pList->detectors[pList->ndet].nIso])); 
			
		    pList->detectors[pList->ndet].nIso++;
		    ret = fgets(temp,200,detsINI);		   
	    }
	    ret = fgets(temp,200,detsINI);	
	    ret = fgets(temp,200,detsINI);		    
	    
	    std::string comma;
        int isoj=0;
        while( temp[0]!='_' && isoj < pList->detectors[pList->ndet].nIso)
        {
            int i = 0;
	        std::istringstream ionizations(temp);

	        while( std::getline(ionizations,comma,',') )
	            pList->detectors[pList->ndet].ionization[isoj][i++] = atof(comma.c_str());

	        isoj++;
   	        ret = fgets(temp,200,detsINI);	
        }   
        if(isoj > 1 && isoj != pList->detectors[pList->ndet].nIso)
        {
	        std::cout << isoj << "Ionizations not listed for each isotope, aborting\n";
	        return 1;
        }
        else if(isoj==1)
        {
            std::cout << "Ionizations will be duplicated for each isotope\n";
            
            for(int j=1; j<pList->detectors[pList->ndet].nIso; j++)
                 pList->detectors[pList->ndet].ionization[j] = pList->detectors[pList->ndet].ionization[0];
        }
        
	    //finished reading in det data		 
		fclose(detsINI);
		
		//only need to calculate background and SM signal once, store in a table for interpolation, stored as events/kg/day/keV
		//get values of bg at relevant energies
		std::cout << "Initializing rates\n"; 		
		rateInit( pList, pList->ndet, -1, &detBackground, pList->detectors[pList->ndet].background);
        std::cout << ".";
    
        //must initialize each flux in turn
        for(int fluxj=0; fluxj< pList->source.numFlux; fluxj++)
        {
            pList->source.nuFluxNorm[fluxj] = 1;
            
		    rateInit( pList, pList->ndet, fluxj,  &SMrate,  pList->detectors[pList->ndet].signalSM[fluxj]);
		    std::cout << ".";
		    
		    //initialization is a little more complicated for the BSM case because of interference terms 
		    if(pList->BSM!=0)
		    {   
		        pList->SMinterference1=1;  pList->SMinterference2=0;
		        rateInit( pList, pList->ndet, fluxj, &BSMrate, pList->detectors[pList->ndet].signalBSM1[fluxj]);
		        
		        if(pList->BSM==3 || pList->BSM==4)
		        {
		                pList->SMinterference2=1; pList->SMinterference1=0;
		                rateInit( pList, pList->ndet, fluxj, &BSMrate,  pList->detectors[pList->ndet].signalBSM2[fluxj]);
		        }
		        pList->SMinterference1=pList->SMinterference2=1;
            }
            else
            {
                pList->SMinterference1=0;  pList->SMinterference2=0;
		        rateInit( pList, pList->ndet, fluxj, &BSMrate, pList->detectors[pList->ndet].signalBSM1[fluxj]);
            }
            std::cout << ".";
        }
        
		std::cout << "\ndone." << std::endl; 
		
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
