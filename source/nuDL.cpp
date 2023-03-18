#include <iostream>
#include <cstring>
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif
#include "fileIO.h"
#include "interactiveInput.h"
#include "SMrate.h"
#include "BSMrate.h"
#include "discLimit.h"
#include "exclusionLimit.h"
#include "calcRates.h"
#include "confInterval.h"
#include "sterileOscillation.h"
double logPoisson(double obs, double expect);
void MCtestDisc(paramList *pList);

int main(int argc, char *argv[])
{
    
    //argument validation and file input
    char filename[100];
    int err=0;
    int mode=0;
    paramList pList;
    
    if(argc==1)
    {
        std::cout << "reading from default config.dat file" << std::endl; 
        sprintf(filename,"config.dat"); 
        mode = readConfigFile(&pList,filename);
    }
    else if(argc>1)
    {
        for(int i=1; i<argc; i++ )
        {
            
            if ( std::strcmp( argv[i], "-c") == 0)
            {
               sprintf(filename,"%s",argv[i+1]);
               mode = readConfigFile(&pList,filename);
               i++;
            }
            else if ( std::strcmp(argv[i], "-i") == 0 )
            {
                std::cout << "Running nuDL in interactive mode" << std::endl;
                mode = interactiveInput(&pList);
            }
            else
            {
                std::cerr << "nuDL: Invalid input" << std::endl << std::endl;
                std::cerr << "Usage: ./nuDL " << std::endl;
                std::cerr << "       default (no flags):" << std::endl; 
                std::cerr << "               (nuDL runs with the default \"config.dat\" parameter file)" << std::endl;
                std::cerr << "       optional flags: " << std::endl; 
                std::cerr << "               -i          (nuDL starts in interactive mode)" << std::endl;
                std::cerr << "               -c file.dat (nuDL runs with the specified parameter file)" << std::endl << std::endl;
                break;
            }
        }
    } 
    
	//pList.printPars();
    
	if ( mode < 1 ) 
    {
        std::cerr << "nuDL: Problem with configuration, aborting" << std::endl;
        return 0;
    }
	
	//print-out rate mode
	if(mode == 1)
    {
        err = calcRates(&pList);
        return 0;
    }   
    
    //discovery limit evolution mode
    if ( mode == 2 )
    {
        discLimitEvolution(&pList, 0);
        return 0;
    }
    
    //discovery limit as a function of mediator mass
    if ( mode == 3 )
    {
        discLimitVsMmed(&pList, 0);
        return 0;
    }
    
    //exclusion limit as a function of mediator mass
    if ( mode == 4 )
    {
        exclusionLimit(&pList, 0);
        return 0;
    }
    
    //print-out total rate above threshold mode
	if(mode == 5)
    {
        err = calcRatesThreshold(&pList);
        return 0;
    }  
    
    //find confidence interval
    if ( mode == 6 )
    {
        //MCtestDisc(&pList);
        confIntVsExposure(&pList);
        return 0;
    }
    
    //discovery limit as a function of threshold
    if ( mode == 7 )
    {
        discLimitVsThresh(&pList, 0);
        return 0;
    }
    
    //standard model significance 
    if ( mode == 8 )
    {
        SMsignificanceExp(&pList);
        return 0;
    }

    //sterile sensitivity mode
    if ( mode == 9 )
    {
        sterileOscillation(&pList);
        return 0;
    }
    
    std::cerr << "nuDL: choose a valid mode" << std::endl;
    return 0;
}
