#include <iostream>
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif
#include "fileIO.h"
#include "CNNSrate.h"
#include "discLimit.h"

int main(int argc, char *argv[])
{

    //argument validation
    char filename[100];
    switch(argc)
    {
        case 1:
        {
            sprintf(filename,"config.dat"); 
            break;
        }
        case 2:
        {
            sprintf(filename,"%s",argv[1]);
            break;
        }
        default:
        {
            std::cout << "Usage: ./Conudl {optional: parameter file name}";
            return 0;
        }
    }		  
    
    int err;
    paramList pList;
    
	//file input
	int mode = readConfigFile(&pList,filename);
	pList.printPars();
	
    if(err) { std::cout << "Flux not normalized, exiting" << std::endl; return 1; }
    
	if ( mode < 1 ) 
    {
        std::cout << "Problem with configuration, aborting" << std::endl;
        return 0;
    }
	
	//print rate mode
	if(mode == 1)
    {
        std::cout << "Er (eV)  " << "SM dN/dE     BSM dN/dE (events/kg/day/keV)" << std::endl;
        for (int i=50; i<1001; i+=10)
            std::cout << (double)i/1e3 << "     " << diffSMrate( (double) i/1e3, &pList, 0) << "      " << diffBSMrate( (double) i/1e3, &pList, 0) << std::endl; 
        
        //std::cout << "total rate: " << pList.rateFunc( pList.detectors[0].ErL, pList.detectors[0].ErU, &pList, 0) << " events/kg/day" << std::endl;
        
        return 0;
    }   
    
    //discovery limit evolution mode
    if ( mode == 2 )
    {
        std::cout << "Starting disc. evolution calculations..." << std::endl;

        discLimit(&pList, 0);

        return 0;
    }
    
    return 0;
}
