//  #include "multinest.h"
#include <iostream>
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif
#include "CNNSrate.h"
#include "detectors.h"
#include "nuFlux.h"

int main(void)
{

    //argument validation
    
    int err;
    
    //file input
    
    //code 
    
    //test
    paramList pList;
    
    //create a detector
    detector geDet;
    char name[] = "#GERMANIUM";
    newDetector(&geDet, name, 1, 0);
    pList.detectors[0] = geDet; pList.ndet++;
    // pList.printPars();
    
    //initialize flux
    err = initFlux(&pList);
    if(err) { std::cout << "Problem with flux, exiting" << std::endl; return 1; }
    
    //set parameters
    pList.qA = -0.5;
    pList.qV = 1;
    
    for (int i=1; i<1001; i+=10)
    {
        std::cout << (double)i << "   " << (double) SMrate((double) i/1e3, &pList, 0) << std::endl; 
    }

/*
//MultiNest starts here
        int pWrap[] = {0,0,0,0};						      // which parameters to have periodic boundary conditions?
        int seed = -1;			   					          // random no. generator seed, if < 0 then take the seed from system clock
      
        int ndims = SET_THIS;       // any combination of mass, sigmaSI, sigmaSIvec, sigmaSD, delta, fn/fp, bn/bp, an/bp and rho_DM, v0, vesc 
        int npar = SET_THIS;                // npar can be greater than ndim if you want to get other values from the loglike function output to file
        double logZero = -DBL_MAX;							  // points with loglike < logZero will be ignored by MultiNest
        int initMPI = 0;								      // initialize MPI routines?, relevant only if compiling with MPI
                                                              // set it to 0 if you want your main program to handle MPI initialization
        int outfile = 1;								      // write output files?
        int updateInt = 100000;								  // update interval (for calling dumper, which isn't used here)
        
    
    nested::run((bool)par.sampling[0],(bool)par.sampling[1],(int)par.sampling[2], par.sampling[6], par.sampling[5],ndims, npar, ndims,       100, updateInt, par.sampling[7],  par.root, seed, pWrap, (bool)par.sampling[3], (bool)par.sampling[4], (bool)outfile, (bool)initMPI, logZero, LogLikedN, dumper, pointer);
       */
       
       //clean up and end
       
       return 0;
}
