//definition of struct which contains source info including fluxes
#include <gsl/gsl_interp.h>
#define SOURCESTRUCT_H

struct sourceStruct {

    char name[50];
    int numFlux;
    double distance;
    
    int isLine[10];
    double lineE[10];
    int isSolar[10];
    double survProb[10];
    double survProbUn[10];
    
    double nuFlux[10];
    double nuFluxUn[10];
    double nuFluxNorm[10];
    
    double *flux_E[10];
    double *flux_N[10];
    int flux_points[10];
    double EnuMax[10];
    
    gsl_interp *nuFluxInterp[10];
	gsl_interp_accel *nuFluxAccel[10];
	
	sourceStruct()
	{
	    distance=1;
	    numFlux=0;
	    
	    for(int i=0;i<10;i++)
	    {
	        nuFlux[i]=0;
	        nuFluxUn[i]=1e-99;
	        nuFluxNorm[i]=0;
	        isSolar[i]=0;
	        isLine[i]=0;
	        survProb[i]=1;
	        survProbUn[i]=1e-99;
	        flux_points[i]=0;
	    }
	}
	
};
