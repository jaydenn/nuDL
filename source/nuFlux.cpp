#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>
#include <string>
#include <fstream>
#include <sstream>
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif  
#ifndef SOURCESTRUCT_H
    #include "sourceStruct.h"
#endif  

#include "physicalConstants.h"

double FIRSTEVALCNS = 0;

gsl_integration_workspace * W;

int nuFluxInit(paramList *pL, std::string sourceName)
{
    //read in source file
    std::ifstream sourceFile("sources.ini");
    
    std::string line;
    std::getline(sourceFile, line);
    
    while( line.compare(sourceName) != 0)
    {
        if(sourceFile.eof())
        {    
            std::cout << "source " << sourceName << " not found\n";
            return -1; 
        } 
        else
            std::getline(sourceFile, line);
    }
    
    //read in component fluxes
    std::string fluxFile;
    int fluxj=0;
    double fluxE,fluxN,lineEnergy;
    std::string plusMinus = "+/-";
    std::string flavor;
    std::getline(sourceFile, line);
    
    while( line[0] != '-' && fluxj < 10 && !sourceFile.eof())
    {
        
        std::istringstream lineStream(line);
        //is this a lineFlux (in energy)?        
        if(line.compare(0,4,"line")==0)
        {
            
            if(!( lineStream >> fluxFile >> lineEnergy >> flavor >> pL->source.nuFlux[fluxj] >> plusMinus >> pL->source.nuFluxUn[fluxj] ))
            {
                std::cout << "error parsing source data (line)\n";
                return -1;
            }
            pL->source.nuFluxUn[fluxj] /= pL->source.nuFlux[fluxj]; //want fractional uncertainty
            pL->source.nuFlux[fluxj] *= pow(HBARC,2);               //convert units
            pL->source.lineE[fluxj] = lineEnergy/1000;
            pL->source.isLine[fluxj] = 1;
        }
        else
        {

            if(!( lineStream >> fluxFile >> flavor >> pL->source.nuFlux[fluxj] >> plusMinus >> pL->source.nuFluxUn[fluxj] ))
            {
                std::cout << " error parsing source data\n";
                return -1;
            }
            
            pL->source.nuFluxUn[fluxj] /= pL->source.nuFlux[fluxj];  //want fractional uncertainty
            pL->source.nuFlux[fluxj] /= pow(pL->source.distance,2);   
            pL->source.nuFlux[fluxj] *= pow(HBARC,2);           //convert units from /cm^2
            
            if( !( lineStream >> pL->source.survProb[fluxj] >> plusMinus >> pL->source.survProbUn[fluxj] ) )
            {
                pL->source.survProb[fluxj]=1;
                pL->source.isSolar[fluxj]=0;
            }
            else
            {
                pL->source.isSolar[fluxj]=1;
                pL->source.nuFluxUn[fluxj] = sqrt( pow(pL->source.nuFluxUn[fluxj],2) + pow(pL->source.survProbUn[fluxj]/pL->source.survProb[fluxj],2) );
            }   
            
            fluxFile.insert (0,  "source/fluxes/");
            std::cout<< "reading " << fluxFile << std::endl;
            std::ifstream flux(fluxFile.c_str());
            pL->source.flux_E[fluxj] = new double [1000];
            pL->source.flux_N[fluxj] = new double [1000];
            
            int i=0;
            int maxPoints=1000;
            while( flux >> fluxE >> fluxN )
            {
                pL->source.flux_E[fluxj][i  ] = fluxE/1000;   //convert to GeV
                pL->source.flux_N[fluxj][i++] = fluxN*1000;
                pL->source.flux_points[fluxj]++;     
                
                if(i==maxPoints)
                {
                    maxPoints+=100;
                    pL->source.flux_E[fluxj] = (double *) realloc(pL->source.flux_E[fluxj], maxPoints * sizeof(double));
                    pL->source.flux_N[fluxj] = (double *) realloc(pL->source.flux_N[fluxj], maxPoints * sizeof(double));
                }       
            }
            flux.close(); flux.clear();
            
            //initialize gsl interpolator for flux
            pL->source.nuFluxInterp[fluxj] = gsl_interp_alloc(gsl_interp_linear, pL->source.flux_points[fluxj]);
            gsl_interp_init(pL->source.nuFluxInterp[fluxj], pL->source.flux_E[fluxj], pL->source.flux_N[fluxj], pL->source.flux_points[fluxj]);
            pL->source.nuFluxAccel[fluxj] = gsl_interp_accel_alloc();
            pL->source.EnuMax[fluxj] = pL->source.flux_E[fluxj][pL->source.flux_points[fluxj]-1];
            
            double norm = gsl_interp_eval_integ (pL->source.nuFluxInterp[fluxj], pL->source.flux_E[fluxj], pL->source.flux_N[fluxj], pL->source.flux_E[fluxj][0], pL->source.flux_E[fluxj][pL->source.flux_points[fluxj]-1], pL->source.nuFluxAccel[fluxj]);
    
            if ( fabs(norm-1.00) > .01)
            {    
                std::cout << "ERROR: flux data in " << fluxFile << " is not properly normalized N = " << norm << std::endl;
                return -1;
            }
        }
        
        if( flavor == "e" )
            pL->source.nuFluxFlav[fluxj] = 1;
        else if( flavor == "ebar")
            pL->source.nuFluxFlav[fluxj] = -1;
        else if( flavor == "mu" )
            pL->source.nuFluxFlav[fluxj] = 2;
        else if( flavor == "mubar" )
            pL->source.nuFluxFlav[fluxj] = -2;
        else if( flavor == "tau" )
            pL->source.nuFluxFlav[fluxj] = 3;
        else if( flavor == "taubar" )
            pL->source.nuFluxFlav[fluxj] = -3;
        else
        {
            std::cout << "flavor of neutrino flux not recognized\n";
            return -1;
        }
        
        fluxj++;
        std::getline(sourceFile, line);

    }
    if(fluxj==10)
        std::cout << "WARNING: max number of flux elements reached (10), ignoring any further components\n";
    
    pL->source.numFlux = fluxj;
    
    sourceFile.close();

    return 0;
    
}

//returns diffNuFlux in GeV per sec, at the point Enu(GeV) ( /cm^2/s/GeV * hc^2)
double nuFlux(double EnuGeV, paramList *pL, int fluxj)
{
    if(EnuGeV < pL->source.EnuMax[fluxj] && EnuGeV > pL->source.flux_E[fluxj][0])
        return pL->source.nuFlux[fluxj] * gsl_interp_eval(pL->source.nuFluxInterp[fluxj], pL->source.flux_E[fluxj], pL->source.flux_N[fluxj], EnuGeV, pL->source.nuFluxAccel[fluxj]);
    else
        return 1e-299;
}

//units GeV/s
double EnuIntegrand0(double EnuGeV, void *pars)
{
	paramList *pList = (paramList*)pars;
	return nuFlux(EnuGeV, pList, pList->fluxj);
}

//units 1/s
double EnuIntegrand1(double EnuGeV, void *pars)
{
	paramList *pList = (paramList*)pars;
	return nuFlux(EnuGeV, pList, pList->fluxj) / EnuGeV ;
}

//units 1/GeV/s
double EnuIntegrand2(double EnuGeV, void *pars)
{
	paramList *pList = (paramList*)pars;
	return nuFlux(EnuGeV, pList, pList->fluxj) / ( EnuGeV*EnuGeV ) ;
}

double fluxIntegral(double ErGeV, paramList *pList, double Mt, int EnuPow, int fluxj, double dist)
{
	int limit = 3000;
	double integral,absErr,tol;
	
	if(FIRSTEVALCNS==0)
	{
		W = gsl_integration_workspace_alloc (3000);
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
		tol=1e-21;
	}
	else if(EnuPow==-2)
	{
		pList->F.function = &EnuIntegrand2;
		tol=1e-20;
	}
			 
	pList->F.params = pList; //yeah, that's not weird..
	
	double EnuMinGeV = 0.5 * (ErGeV + sqrt( pow(ErGeV,2) + 2*ErGeV*Mt ) );
    
    pList->fluxj = fluxj;
    if(pList->source.isLine[fluxj]==1)
    {
        if(EnuMinGeV < pList->source.lineE[fluxj] )
            integral = pList->source.nuFlux[fluxj] * pow( pList->source.lineE[fluxj], EnuPow);
        else
            integral = 1e-199;
    }
    else
        gsl_integration_qag(&(pList->F), EnuMinGeV, pList->source.EnuMax[fluxj], tol, 5e-4, limit, 2, W, &integral, &absErr); 
    
    if (integral < 0)        
        return 1e-299;
    else    
    	return pList->source.nuFluxNorm[fluxj]*integral / pow(dist,2);
		
}


//units GeV/s
double EnuIntegrandOsc0(double EnuGeV, void *pars)
{
	paramList *pList = (paramList*)pars;
	return (1.0-pList->ss2Theta14*pow(sin(pList->delMsqGeV*pList->dist/GEVMETER/4.0/EnuGeV),2))*nuFlux(EnuGeV, pList, pList->fluxj);
}

//units 1/s
double EnuIntegrandOsc1(double EnuGeV, void *pars)
{
	paramList *pList = (paramList*)pars;
	return (1.0-pList->ss2Theta14*pow(sin(pList->delMsqGeV*pList->dist/GEVMETER/4.0/EnuGeV),2))*nuFlux(EnuGeV, pList, pList->fluxj) / EnuGeV ;
}

//units 1/GeV/s
double EnuIntegrandOsc2(double EnuGeV, void *pars)
{
	paramList *pList = (paramList*)pars;
	return (1-pList->ss2Theta14*pow(sin(pList->delMsqGeV*pList->dist/GEVMETER/4.0/EnuGeV),2))*nuFlux(EnuGeV, pList, pList->fluxj) / ( EnuGeV*EnuGeV ) ;
}

//se
double fluxIntegralOsc(double ErGeV,  paramList *pList, double Mt, int EnuPow, int fluxj, double dist)
{
	int limit = 3000;
	double integral,absErr,tol;

	if(FIRSTEVALCNS==0)
	{
		W = gsl_integration_workspace_alloc (3000);
		FIRSTEVALCNS=1;
	}
	
	if(EnuPow==0)
	{
		pList->F.function = &EnuIntegrandOsc0;
		tol=1e-21;
	}
	else if(EnuPow==-1)
	{
		pList->F.function = &EnuIntegrandOsc1;
		tol=1e-19;
	}
	else if(EnuPow==-2)
	{
		pList->F.function = &EnuIntegrandOsc2;
		tol=1e-18;
	}
			 
	pList->F.params = pList; //yeah, that's not weird..
    pList->dist = dist;
    
	double EnuMinGeV = 0.5 * (ErGeV + sqrt( pow(ErGeV,2) + 2*ErGeV*Mt ) );
    
    pList->fluxj = fluxj;
        
    if(pList->source.isLine[fluxj]==1)
    {
        if(EnuMinGeV < pList->source.lineE[fluxj] )
            integral = pList->source.nuFlux[fluxj] * pow( pList->source.lineE[fluxj], EnuPow);
        else
            integral = 1e-199;
    }
    else
        gsl_integration_qag(&(pList->F), EnuMinGeV, pList->source.EnuMax[fluxj], tol, 5e-4, limit, 2, W, &integral, &absErr); 

    if (integral < 0)        
        return 1e-299;
    else    
    	return pList->source.nuFluxNorm[fluxj] * integral / pow(dist,2);
		
}

