#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif  

//data from Joel's notebook
double flux_E[69] = { 0., 0.00001, 0.00002, 0.000035, 0.00004, 0.00007, 0.0001, 0.00013, 0.00016, 0.000165, 0.00018, 0.000215, 0.00023, 0.00028, 0.00033, 0.000335, 0.00035, 0.00039, 0.0004, 0.000435, 0.00044, 0.0005, 0.0007, 0.0009, 0.001, 0.001185, 0.00119, 0.00125, 0.0013, 0.0015, 0.0017, 0.0018, 0.0019, 0.002, 0.00225, 0.0025, 0.00275, 0.003, 0.00325, 0.0035, 0.00375, 0.004, 0.00425, 0.0045, 0.00475, 0.005, 0.00525, 0.0055, 0.00575, 0.006, 0.00625, 0.0065, 0.00675, 0.007, 0.00725, 0.0075, 0.00775, 0.008, 0.0081, 0.0082, 0.00825, 0.0083, 0.0084, 0.0085, 0.00875, 0.009, 0.00925, 0.0095, 0.01 };
double flux_N[69] = {2.95612, 11.496, 44.7068, 125.951, 52.5788, 146.894, 268.835, 406.966, 571.831, 484.2, 555.493, 692.138, 605.992, 747.094, 886.71, 638.669, 605.992, 656.492, 600.051, 652.036, 421.818, 450.039, 466.377, 479.744, 461.921, 415.877, 332.702, 305.967, 262.894, 236.159, 219.821, 210.909, 201.997, 193.086, 160.41, 133.675, 113.029, 94.6121, 79.6108, 64.9065, 52.2817, 42.0333, 33.1216, 25.5467, 19.6056, 15.5954, 12.1941, 9.16415, 7.15903, 5.49552, 4.01024, 3.01511, 2.22791, 1.55954, 0.992164, 0.637183, 0.399539, 0.201997, 0.13219, 0.0745608, 0.0613419, 0.0481229, 0.0414392, 0.035201, 0.0191601, 0.00831754, 0.00326761, 0.00207939, 0.};
const int flux_points = 69;

//double flux_E[] = {.000,.00016845,.00024648,.00035141,.00051864,.0008492,.00116954,.001503,.00208795,.00297672,.00428068,.00584437,.00757548,.00892807,.0105219,.01167139,.01272424,.0136333,.014480,.01498602,.01564164,.01604703,.01631117,.01631598,.01658544};
//double flux_N[] = {1.05911035,3.11556519,5.9498291,11.36418555,21.70120703,54.27990234,98.23500781,135.7049375,245.58442188,377.794125,581.1490625,720.232625,681.2650625,548.1728125,375.002625,195.808375,91.76617188,27.90096289,4.19898828,0.82834375,0.04970127,0.00459756,0.00042531,0.00000293,0.00000011};
//const int flux_points = 25;

double FIRSTEVALCNS = 0;
gsl_integration_workspace * W;

int initFlux(paramList *pList)
{
    pList->nuFluxInterp = gsl_interp_alloc(gsl_interp_linear, flux_points);
    gsl_interp_init(pList->nuFluxInterp, flux_E, flux_N, flux_points);
    pList->nuFluxAccel = gsl_interp_accel_alloc();

    pList->EnuMax = flux_E[flux_points-1];
    
    double norm = gsl_interp_eval_integ (pList->nuFluxInterp, flux_E, flux_N, flux_E[0], flux_E[flux_points-1], pList->nuFluxAccel);
    
    if ( fabs(norm-1.00) > .01)
        return 1;
    else
        return 0;
}

//returns diffNuFlux in GeV per sec, at the point Enu(GeV) ( /cm^2/s/GeV * hc^2)
double nuFlux(double EnuGeV, paramList *pList)
{
    if(EnuGeV < pList->EnuMax && EnuGeV > flux_E[0])
        return pList->nuFlux * gsl_interp_eval(pList->nuFluxInterp, flux_E, flux_N, EnuGeV, pList->nuFluxAccel);
    else
        return 1e-99;
}


//units GeV/s
double EnuIntegrand0(double EnuGeV, void *pars)
{
	paramList *pList = (paramList*)pars;
	return nuFlux(EnuGeV, pList);
}

//units 1/s
double EnuIntegrand1(double EnuGeV, void *pars)
{
	paramList *pList = (paramList*)pars;
	return nuFlux(EnuGeV, pList) /	EnuGeV ;
}

//units 1/GeV/s
double EnuIntegrand2(double EnuGeV, void *pars)
{
	paramList *pList = (paramList*)pars;
	return nuFlux(EnuGeV, pList) / ( EnuGeV*EnuGeV ) ;
}

double fluxIntegral(double ErGeV,  paramList *pList, double Mt, int EnuPow)
{
	int key = 2;
	int limit = 2000;
	double integral,absErr,tol;
	
	if(FIRSTEVALCNS==0)
	{
		W = gsl_integration_workspace_alloc (2000);
		FIRSTEVALCNS=1;
	}
	
	if(EnuPow==0)
	{
		pList->F.function = &EnuIntegrand0;
		tol=1e-22;
	}
	else if(EnuPow==-1)
	{
		pList->F.function = &EnuIntegrand1;
		tol=1e-20;
	}
	else if(EnuPow==-2)
	{
		pList->F.function = &EnuIntegrand2;
		tol=1e-18;
	}
			 
	pList->F.params = pList; //yeah, that's not weird..
	
	double EnuMinGeV = sqrt( ErGeV * Mt / 2);

	gsl_integration_qag(&(pList->F), EnuMinGeV, pList->EnuMax, tol, 1e-3, limit, 2, W, &integral, &absErr); 
    //std::cout << ErGeV  << " " << EnuIntegrand2(ErGeV,(void *)pList) << " " << EnuIntegrand1(ErGeV,(void *)pList) << " " << EnuIntegrand0(ErGeV,(void *)pList) << std::endl;
	
	return integral;	
}

