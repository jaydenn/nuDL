#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_multimin.h>
#include <string.h>
#include <gsl/gsl_cdf.h>
#ifndef DETECTORSTRUCT_H
	#include "detectorStruct.h"
#endif
#ifndef PARAMETERSTRUCT_H
	#include "parameterStruct.h"
#endif	
#ifndef DETECTORFUNCTIONS_H
	#include "detectorFunctions.h"
#endif
#include "likelihood.h"
#include "monteCarlo.h"
#include "nuRate.h"
#include "SMrate.h"
#include "BSMrate.h"


double L(const gsl_vector *v, void *params)
{

    paramList *pL = (paramList *)params;
    double l = 0;
    
    pL->signalNorm = gsl_vector_get(v, 0);
    pL->detectors[pL->detj].BgNorm = fabs(gsl_vector_get(v, 1));
    l += 0.5/pow(pL->detectors[pL->detj].BgUn,2) * pow(pL->detectors[pL->detj].BgNorm - 1,2);
    
    for(int i=0; i < pL->source.numFlux; i++)
    {
        pL->source.nuFluxNorm[i] = fabs(gsl_vector_get(v, i+2));
        l += 0.5/pow(pL->source.nuFluxUn[i],2) * pow(pL->source.nuFluxNorm[i] - 1,2);
    }
    
    l -= logLikelihoodSM(pL);
}

double Lmu(const gsl_vector *v, void *params)
{

    paramList *pL = (paramList *)params;
    double l = 0;

    pL->detectors[pL->detj].BgNorm  = gsl_vector_get(v, 0);
    l += 0.5/pow(pL->detectors[pL->detj].BgUn,2) * pow(pL->detectors[pL->detj].BgNorm - 1,2);

    std::cout <<  "\n 98t \n" ;
    for(int i=0; i < pL->source.numFlux; i++)
    {
        std::cout <<  "\n y " << i;
        pL->source.nuFluxNorm[i] = fabs(gsl_vector_get(v, i+1));
        l += 0.5/pow(pL->source.nuFluxUn[i],2) * pow(pL->source.nuFluxNorm[i] - 1,2);
    }
    
    l -= logLikelihoodSM(pL);
    
    return l;

}

double findMaxL(paramList *pL)
{
    size_t iter = 0;
    int status;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_multimin_function my_func;
    
    my_func.f = L;
    my_func.params = (void *)pL;

    //start point and step size
    gsl_vector *x,*dx;
    my_func.n = 2 + pL->source.numFlux;
    
    x = gsl_vector_alloc ( my_func.n );
    dx = gsl_vector_alloc ( my_func.n );
    gsl_vector_set (x, 0, pL->signalNorm);
    gsl_vector_set(dx, 0, pL->signalNorm/10);

    for(int i=1; i < my_func.n; i++)
    {
        gsl_vector_set (x, i, 1.0);
        gsl_vector_set(dx, i, .05);
    }
    
    std::cout <<  " 6t " <<     pL->signalNorm << std::endl;
    gsl_multimin_fminimizer_set (s, &my_func, x, dx);
    std::cout <<  " 7t " <<     pL->signalNorm << std::endl;
    
    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate (s);       
    //        cout << "   " <<iter << " " <<  gsl_vector_get (s->x, 0) << " " <<  gsl_vector_get (s->x, 1) << " " <<  gsl_vector_get (s->x, 2)*3e5 << " " <<  gsl_vector_get (s->x, 3)*3e5 << " " << s->fval << " " << gsl_multimin_fminimizer_size(s) << endl; 
    }
    while (iter < 900 && gsl_multimin_fminimizer_size(s)>.0001); //s->fval>1e-2);

   if(iter==900)
        std::cout << "non-convergence size = " << gsl_multimin_fminimizer_size(s)/s->fval << " > .0001  " << std::endl;
    
    double LS =  s->fval;
    pL->signalNorm = gsl_vector_get(s->x, 1);
    
    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free (dx);
    
    return LS;
}


double findMaxLMu(paramList *pL)
{
    size_t iter = 0;
    int status;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_multimin_function my_func;

    
    my_func.f = Lmu;
    my_func.params = (void *)pL;

    //start point
     gsl_vector *x,*dx;
    my_func.n = 1 + pL->source.numFlux;
    
    x = gsl_vector_alloc (my_func.n);
    dx = gsl_vector_alloc (my_func.n);
    for(int i=0; i < my_func.n; i++)
    {
        gsl_vector_set (x, i, 1.0);
        gsl_vector_set(dx, i, .05);
    }
    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, my_func.n);

    gsl_multimin_fminimizer_set (s, &my_func, x, dx);
    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate (s);
       //std::cout << "       " << iter << "  c " << " " <<  gsl_vector_get (s->x, 0) << " " <<  gsl_vector_get (s->x, 1) << " " <<  gsl_vector_get (s->x, 2) << " " << s->fval << " " << gsl_multimin_fminimizer_size(s) << std::endl; 
    }
    while (iter < 900 && gsl_multimin_fminimizer_size(s)>.0001);///s->fval > 1e-2 && !status);
    if(iter==900)
        std::cout << "non-convergence size = " << gsl_multimin_fminimizer_size(s)/s->fval  << " > .0001  " <<  std::endl;
    
    double Lmu = s->fval;
    std::cout <<  " 6t " <<     pL->signalNorm << std::endl;
    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free (dx);
            
    return Lmu;
}

//modified test statistic for mu >= 0, returns log of test stat
double tMu(paramList *pL)
{	 

    double maxLmu = findMaxLMu( pL );
    std::cout << "maxLmu = " << maxLmu << std::endl;
    double maxL = findMaxL( pL );
    
    if( pL->signalNorm < 0)
    {   
        pL->signalNorm = 0;
        maxL = findMaxLMu( pL );
    }
    
    std::cout << maxLmu << " " << maxL << std::endl;
    
    
    return ( maxL - maxLmu );  //no -ve because functions return -loglike

}

double pValue(double tMu)
{

    double cdf = 2 * gsl_cdf_ugaussian_P(sqrt(tMu)) -1;
    return 1 - cdf;
}

double my_tMu(const gsl_vector *v, void *params)
{

    double simSeed = 0;
    paramList *pL = (paramList *)params;
    pL->signalNorm  = gsl_vector_get(v, 0);
    std::cout <<  " t " <<     pL->signalNorm << std::endl;
    
    double pVal = pValue( tMu(pL) );
    
    return pow(pVal - 0.9, 2);  //90% confidence exclusion
    
}


double findtMu90(paramList *pL)
{

    size_t iter = 0;
    int status;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;

    gsl_vector *x,*dx;
    gsl_multimin_function my_func;

    my_func.n = 1;
    my_func.f = my_tMu;
    my_func.params = (void *)pL;

    //start point
    x = gsl_vector_alloc (1);
    dx = gsl_vector_alloc (1);
    
    gsl_vector_set (x, 0, pL->signalNorm);
    gsl_vector_set(dx, 0, pL->signalNorm/10);

    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, 1);

    gsl_multimin_fminimizer_set (s, &my_func, x, dx);

    do
    {
        status = gsl_multimin_fminimizer_iterate (s);
        iter++;
      //  std::cout << iter << " c = " << gsl_vector_get(s->x,0) << ",  f = " << s->fval << ", q= " << 0.5*pow(sqrt( s->fval)+1.282,2)+1 << "  size = " << gsl_multimin_fminimizer_size(s) << std::endl;
    }
    while (iter < 1200 && s->fval > .0002 && !status); 
    if(iter==1200)
        std::cout << "non-convergence (sig-1.28)^2 = " << s->fval << " > .0002" << std::endl;
        
    double c = gsl_vector_get(s->x, 0);
    
    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free (dx);

    return c;
}

void confIntervalSM(paramList *pL)
{
    int detj = 0;
    //add
    pL->detj = detj;    

    //find a good starting point
    double BG  = intBgRate(pL->detectors[detj], pL->detectors[detj].ErL, pL->detectors[detj].ErU);         
    double SM  = intSMrate( pL->detectors[detj].ErL, pL->detectors[detj].ErU, pL, detj);

    double mu = pL->signalNorm = BG/SM/5;
    std::cout << "starting mu = " << mu << std::endl;

    pL->BSM = 0;
    generateBinnedData( pL, pL->detj, 1, 0);
    
    double muUp = findtMu90(pL);

    //print out result
    std::cout << mu << std::endl;
            
    pL->signalNorm = mu/10;      //update guess to find other solution
        
    double muDown = findtMu90(pL);

    //print out result
    std::cout << "90\% confidence interval of mu = " << muDown << " - " << muUp << std::endl;
    
}

double confIntVsBlah(paramList *pL)
{
    /*if(pL->elecScat)
        sprintf(filename, "%sconfIntE_%c%c_BSM%d.dat",pL->root,pL->detectors[0].name[0],pL->detectors[0].name[1],pL->BSM);
    else
        sprintf(filename, "%sconfIntN_%c%c_BSM%d.dat",pL->root,pL->detectors[0].name[0],pL->detectors[0].name[1],pL->BSM);
    
    outfile.open(filename,std::ios::out);
    std::cout << "writing output to: " << filename << std::endl;
    */
    return 0.0;
}

