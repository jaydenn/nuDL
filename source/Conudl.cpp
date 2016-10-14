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
                std::cout << "Running Conudl in interactive mode" << std::endl;
                mode = interactiveInput(&pList);
            }
            else
            {
                std::cout << "Invalid input" << std::endl << std::endl;
                std::cout << "Usage: ./Conudl " << std::endl;
                std::cout << "       default (no flags):" << std::endl; 
                std::cout << "               (Conudl runs with the default \"config.dat\" parameter file)" << std::endl;
                std::cout << "       optional flags: " << std::endl; 
                std::cout << "               -i          (Conudl starts in interactive mode)" << std::endl;
                std::cout << "               -c file.dat (Conudl runs with the specified parameter file)" << std::endl << std::endl;
                break;
            }
        }
    } 
    
	//pList.printPars();
    
	if ( mode < 1 ) 
    {
        std::cout << "Problem with configuration, aborting" << std::endl;
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
        confIntervalSM(&pList);
        return 0;
    }
    
    std::cout << "Choose a valid mode" << std::endl;
    return 0;

}
