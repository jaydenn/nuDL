//definition of parameter structs
#include "detectorStruct.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>
#define PARAMETERSTRUCT_H

struct paramList {
    
    //arrays for reactor flux data
    gsl_interp *nuFlux;
    gsl_interp_accel *nuFluxAccel;
    double nuFluxNorm, EnuMax;
    
    //setup for rate integration
    gsl_function F;
     
    double Er;
    double A;
    double qA;
    double qAn;
    double qAp;
    double qV;
    double qVn;
    double qVp;
    
    double ndet;
    detector detectors[10];
    
    void printPars()
    {
        
        printf("Detectors:");
        for(int i=0;i<ndet;i++)
            detectors[i].printDetSpecs();
    }
    
    paramList()
    {
        ndet=0;
        nuFluxNorm=1;
        EnuMax=0;
    }
};
    
