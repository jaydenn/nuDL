#include <iostream>
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif
#include "CNNSrate.h"
#include "discLimit.h"
#include "nuFlux.h"

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
        
    //test
    paramList pList;
    
	//file input
	//readConfigFile(filename, &pList);
	
    //create a detector    
    char name[] = "#GERMANIUM";
	double exposure = 1;
    err = pList.newDetector(name, exposure);
	if(err) { std::cout << "Could not create detector, exiting" << std::endl; return 1; }
	
    // pList.printPars();
    
    //initialize reactor flux
    err = initFlux(&pList);
    if(err) { std::cout << "Problem with flux, exiting" << std::endl; return 1; }
    
    SMrateInit(&pList, 0);
    pList.rateFunc = &intSMrate;    
    discLimit(&pList, 0);
       
    for (int i=50; i<1001; i+=10)
    {
        std::cout << (double)i << "   " << pList.detectors[0].exposure * intCNNSrate( (double) i/1e3, (double) (i+10)/1e3, &pList, 0) << "    " << 100*pList.detectors[0].exposure*1e-2 << std::endl; 
    }
       
    return 0;
}
