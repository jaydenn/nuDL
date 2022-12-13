#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_multimin.h>
#include <string.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sort_vector.h>
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
#include "sterileRate.h"
#include "SMrate.h"
#include <gsl/gsl_errno.h>

double my_LS_sterile(const gsl_vector *v, void *params)
{

    paramList *pL = (paramList *)params;
    double l = 0;

    pL->signalNorm = gsl_vector_get(v, 0);
    pL->ss2Theta14 = gsl_vector_get(v, 0);
    
    pL->detectors[pL->detj].BgNorm = fabs(gsl_vector_get(v, 1));
    l += 0.5/pow(pL->detectors[pL->detj].BgUn,2) * pow(pL->detectors[pL->detj].BgNorm - 1,2);
    
    for(int i=0; i < pL->source.numFlux; i++)
    {
        pL->source.nuFluxNorm[i] = fabs(gsl_vector_get(v, i+2));
        l += 0.5/pow(pL->source.nuFluxUn[i],2) * pow(pL->source.nuFluxNorm[i] - 1,2);
    }
    l -= logLikelihoodSterile(pL);
    //std::cout << " like "  << l << " " << pL->signalNorm << " " << pL->source.nuFluxNorm[0]  << "\n";
    return l;

}

double my_L0_sterile(const gsl_vector *v, void *params)
{

    paramList *pL = (paramList *)params;
    double l = 0;

    pL->detectors[pL->detj].BgNorm  = gsl_vector_get(v, 0);
    l += 0.5/pow(pL->detectors[pL->detj].BgUn,2) * pow(pL->detectors[pL->detj].BgNorm - 1,2);
    
    for(int i=0; i < pL->source.numFlux; i++)
    {
        pL->source.nuFluxNorm[i] = fabs(gsl_vector_get(v, i+1));
        l += 0.5/pow(pL->source.nuFluxUn[i],2) * pow(pL->source.nuFluxNorm[i] - 1,2);
    }
 
    l -= logLikelihoodSM(pL);
    //std::cout << l << std::endl;
    return l;

}

double findMaxLS_sterile(paramList *pL)
{

    size_t iter = 0;
    int status;
    status =0;
    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_multimin_function my_func;
    
    my_func.f = my_LS_sterile;
    my_func.params = (void *)pL;

    //start point and step size
    gsl_vector *x,*dx;
    my_func.n = 2 + pL->source.numFlux;

    x = gsl_vector_alloc ( my_func.n );
    dx = gsl_vector_alloc ( my_func.n );
    gsl_vector_set (x, 0,.5);
    gsl_vector_set(dx, 0,  .05);
    
    for(int i=1; i < my_func.n; i++)
    {
        gsl_vector_set (x, i, 1.0);
        gsl_vector_set(dx, i, .01);
    }
    
    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, my_func.n);

    gsl_multimin_fminimizer_set (s, &my_func, x, dx);

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate (s);
        
        //std::cout << "       " <<iter << " " <<  gsl_vector_get (s->x, 0) << " " <<  gsl_vector_get (s->x, 1) << " " << gsl_vector_get (s->x, 2) << " L " << s->fval << std::endl; 
    }
    while (iter < 1000 && gsl_multimin_fminimizer_size(s)/s->fval > .00001);
    
   if(iter==1000)
        std::cout << "LS non-convergence size = " << gsl_multimin_fminimizer_size(s)/s->fval << " > .0005  " << std::endl;
    
    double LS =  s->fval;
    pL->signalNorm = gsl_vector_get(s->x, 1);
    
    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free (dx);
    
    return -LS;
}


double findMaxL0_sterile(paramList *pL)
{

    size_t iter = 0;
    int status;
    
    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_multimin_function my_func;

    my_func.f = my_L0_sterile;
    my_func.params = (void *)pL;

    //start point
    gsl_vector *x,*dx;
    my_func.n = 1 + pL->source.numFlux;
    
    x = gsl_vector_alloc (my_func.n);
    dx = gsl_vector_alloc (my_func.n);
    
    for(int i=0; i < my_func.n; i++)
    {
        gsl_vector_set (x, i, 1.0);
        gsl_vector_set(dx, i, .01);
    }
    
    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, my_func.n);

    gsl_multimin_fminimizer_set (s, &my_func, x, dx);

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate (s);
        //std::cout << "       " << iter << " " <<  gsl_vector_get (s->x, 0) << " " << gsl_vector_get (s->x, 1) << " " << s->fval << std::endl; 
    }
    while (iter < 800 && gsl_multimin_fminimizer_size(s)/s->fval > 1e-9 && !status);
    if(iter==800)
        std::cout << "L0 non-convergence size = " << gsl_multimin_fminimizer_size(s)/s->fval  << " > 1e-9 " <<  std::endl;

    double L0 = s->fval;
    
    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free (dx);
    
    return -L0;
}

//statistic for discovery
double q0_sterile(paramList *pL)
{	 

    double maxL;
    double maxL0 = findMaxL0_sterile( pL );
    //if using asimov just set parameters to MLE
    if(pL->asimov==1)
    {
        pL->detectors[pL->detj].BgNorm = 1;
        for(int fluxj=0; fluxj < pL->source.numFlux; fluxj++)
            pL->source.nuFluxNorm[fluxj] = 1.0;
        maxL = logLikelihoodSterile(pL);
    }
    else
        maxL = findMaxLS_sterile( pL );
    
    double q = - 2 * ( maxL0 - maxL ); 
    //std::cout << "returning q= " << q << " = -2 *( L0(" << maxL0 << ") - maxL(" << maxL <<") )"<< std::endl;
    if ( pL->signalNorm >= 0 && q > 0 ) //this is to catch roundoff error, but it could hide bugs
    {
        return - 2 * ( maxL0 - maxL );  
    }
    else
    {
        //std::cout << "returning q=0 " << pL->signalNorm << " " << maxL0 << " " << maxL << std::endl;
        return 0;
    } 
}

//searching for the mu which gives a median significance of 4.38sigma
double my_q0_sterile(const gsl_vector *v, void *params)
{

    double simSeed = 0;
    paramList *pL = (paramList *)params;    
    
    for(int i=0; i < pL->source.numFlux; i++)
        pL->source.nuFluxNorm[i] = 1;
        
    pL->detectors[pL->detj].BgNorm = 1;
    pL->ss2Theta14 = gsl_vector_get(v, 0);
    if( pL->ss2Theta14 < 0 )
        return 1e99;
    generateBinnedDataSterile( pL, pL->detj, 1, simSeed);
    
    return pow( sqrt(q0_sterile(pL)) - 4.28, 2);  //arbitrary function with a minima at "3 sigma 90% of the time"
}


double findCoeff3sig_sterile(paramList *pL)
{

    size_t iter = 0;
    int status;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;

    gsl_vector *x,*dx;
    gsl_multimin_function my_func;

    my_func.n = 1;
    my_func.f = my_q0_sterile;
    my_func.params = (void *)pL;

    //start point
    x = gsl_vector_alloc (1);
    dx = gsl_vector_alloc (1);
    
    gsl_vector_set (x, 0, pL->ss2Theta14);
    gsl_vector_set(dx, 0, pL->ss2Theta14);

    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, 1);

    gsl_multimin_fminimizer_set (s, &my_func, x, dx);

    do
    {
        status = gsl_multimin_fminimizer_iterate (s);
        iter++;
        std::cout << iter << " " << gsl_vector_get(s->x,0) << ", q' = " << s->fval << ", size " << gsl_multimin_fminimizer_size(s) << std::endl;
    }
    while (iter < 30 && s->fval > .0014 && !status); //under 1% error in 4.28 sigma value
        
    pL->ss2Theta14 = gsl_vector_get(s->x, 0);
    
    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free (dx);

    if(iter==30)
    {
        double approxError = sqrt(s->fval)/4.28*100;
        std::cout << "WARNING: non-convergence, sigma - 4.28 = " << s->fval << " " << approxError << "% error" << std::endl;
        if (approxError > 2)
            return 0;
        else
            return pL->ss2Theta14;
    }
    else
        return pL->ss2Theta14;
}

void sterileOscillation(paramList *pList)
{
    std::ofstream outfile;
    outfile.open("results/sterile.dat",std::ios::out);
    outfile   << "delMsq (GeV)  ss2Theta14  significance \n";
    
    //dont do this
    gsl_set_error_handler_off();
    
    pList->delMsqGeV = 9e-20;
    double q0;
    //scan mass diff and find theta14 that gives 3sigma
    while (pList->delMsqGeV<1e-16)
    {
        //findCoeff3sig_sterile(pList);
        pList->ss2Theta14 = 0.01;
        while (pList->ss2Theta14<1)
        {
            rateInit( pList, 0, 0, &sterileRate, pList->detectors[0].signalBSM1[0]);
            generateBinnedDataSterile( pList, 0, 1, 1);
            q0 = q0_sterile(pList);

            std::cout << pList->delMsqGeV  << " " << pList->ss2Theta14 << " " << sqrt(q0) <<  std::endl;
            outfile   << pList->delMsqGeV  << " " << pList->ss2Theta14 << " " << sqrt(q0) << "\n";
            
            pList->ss2Theta14 *= 1.28;    
        }
        
        pList->delMsqGeV  *= 1.28;
    }
    outfile.close();
}
