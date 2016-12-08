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

//global likelihood
double Lg(const gsl_vector *v, void *params)
{

    paramList *pL = (paramList *)params;
    double l = 0;
    
    pL->signalNorm = gsl_vector_get(v, 0);
    
    pL->detectors[pL->detj].BgNorm = fabs(gsl_vector_get(v, 1));
    l += -0.5/pow(pL->detectors[pL->detj].BgUn,2) * pow(pL->detectors[pL->detj].BgNorm - 1,2);
    
    for(int i=0; i < pL->source.numFlux; i++)
    {
        pL->source.nuFluxNorm[i] = fabs(gsl_vector_get(v, i+2));
        l += -0.5/pow(pL->source.nuFluxUn[i],2) * pow(pL->source.nuFluxNorm[i] - 1,2);
    }
    
    l += logLikelihoodSM(pL);
    
    return -l;
}

double Lmu(const gsl_vector *v, void *params)
{

    paramList *pL = (paramList *)params;
    double l = 0;

    pL->detectors[pL->detj].BgNorm  = gsl_vector_get(v, 0);
    l += -0.5/pow(pL->detectors[pL->detj].BgUn,2) * pow(pL->detectors[pL->detj].BgNorm - 1,2);

    for(int i=0; i < pL->source.numFlux; i++)
    {
        pL->source.nuFluxNorm[i] = fabs(gsl_vector_get(v, i+1));
        l += -0.5/pow(pL->source.nuFluxUn[i],2) * pow(pL->source.nuFluxNorm[i] - 1,2);
    }
    
    l += logLikelihoodSM(pL);
    
    return -l;

}

double findMaxL(paramList *pL)
{
    size_t iter = 0;
    int status;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_multimin_function my_func;
    
    my_func.f = Lg;
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
    
    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, my_func.n);
    
    gsl_multimin_fminimizer_set (s, &my_func, x, dx);
    
    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate (s);       
        //std::cout << " Ll  " << iter << " " <<  gsl_vector_get (s->x, 0) << " " <<  gsl_vector_get (s->x, 1) << " " << gsl_vector_get (s->x, 2) << " " <<  s->fval << " " << gsl_multimin_fminimizer_size(s) << std::endl; 
    }
    while (iter < 900 && gsl_multimin_fminimizer_size(s)>.0001); //s->fval>1e-2);

   if(iter==900)
        std::cout << "non-convergence size = " << gsl_multimin_fminimizer_size(s)/s->fval << " > .0001  " << std::endl;
    
    double LS =  s->fval;
    pL->signalNorm = gsl_vector_get(s->x, 1);
    
    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free (dx);
    
    return -LS;
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
        gsl_vector_set(dx, i, .01);
    }
    
    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, my_func.n);

    gsl_multimin_fminimizer_set (s, &my_func, x, dx);
    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate (s);
        //std::cout << "   Lmu    " << iter << "  c " << " " <<  gsl_vector_get (s->x, 0) << " " <<  gsl_vector_get (s->x, 1) << " " << s->fval << " " << gsl_multimin_fminimizer_size(s) << std::endl; 
    }
    while (iter < 900 && gsl_multimin_fminimizer_size(s)>.00001);///s->fval > 1e-2 && !status);
    if(iter==900)
        std::cout << "non-convergence size = " << gsl_multimin_fminimizer_size(s)/s->fval  << " > .00001  " <<  std::endl;
    
    double Lmu = s->fval;

    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free (dx);
            
    return -Lmu;
}

//modified test statistic for mu >= 0, returns log of test stat
double tMu(paramList *pL)
{	 

    double maxLmu = findMaxLMu( pL );
    double maxL;
    //if using asimov just set parameters to MLE
    if(pL->asimov)
    {
        pL->signalNorm = 1;
        pL->detectors[pL->detj].BgNorm = 1;
        for(int fluxj=0; fluxj < pL->source.numFlux; fluxj++)
            pL->source.nuFluxNorm[fluxj] = 1.0;
        maxL = logLikelihoodSM(pL);
    }
    else
        maxL = findMaxL( pL );
    
    //use below for modified test statistic, p-val gets more complicated
    /*if( pL->signalNorm < 0)
    {   
        pL->signalNorm = 0;
        maxL = findMaxLMu( pL );
    }*/
    double t = -2*( maxLmu - maxL);
    if (t > 0)
        return t;
    else 
        return 0;

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
    
    double pVal = pValue( tMu(pL) );
    
    return pow(pVal - 0.1, 2);  //90% confidence exclusion
    
}

double q0SM(paramList *pL)
{
    pL->signalNorm = 0;
    double maxL0 = findMaxLMu( pL );
    double maxL;
    //if using asimov just set parameters to MLE
    if(pL->asimov)
    {
        pL->signalNorm = 1;
        pL->detectors[pL->detj].BgNorm = 1;
        for(int fluxj=0; fluxj < pL->source.numFlux; fluxj++)
            pL->source.nuFluxNorm[fluxj] = 1.0;
        maxL = logLikelihoodSM(pL);
    }
    else
        maxL = findMaxL( pL );
    
    if( pL->signalNorm >= 0 )
        return - 2 * ( maxL0 - maxL );  
    else
        return 0;
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
    gsl_vector_set(dx, 0, pL->signalNorm/20);

    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, 1);

    gsl_multimin_fminimizer_set (s, &my_func, x, dx);

    do
    {
        status = gsl_multimin_fminimizer_iterate (s);
        iter++;
       // std::cout << iter << " mu = " << gsl_vector_get(s->x,0) << ",  f = " << s->fval << "  size = " << gsl_multimin_fminimizer_size(s) << std::endl;
    }
    while (iter < 100 && s->fval > 1e-6 && !status); //gives ~1% accuracy in confidence
    if(iter==100)
    {
        std::cout << "non-convergence (p-0.9)^2 = " << s->fval << " > 1e-6" << std::endl;
        return 0;   
    }
        
    double mu = gsl_vector_get(s->x, 0);
    
    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free (dx);

    return mu;
}

double *confIntervalSM(paramList *pL)
{
    pL->detj = 0;    
    double *CI = new double[2];
    int nEvents = generateBinnedData( pL, pL->detj, 1, 0);
    
    double muUp,muDown;
    double pVal = 0;
    double guess = 2 + 10*pL->source.isSolar[0];
    while(guess > 1e-3 &&  pVal < .1)
    {
        guess*=.97;
        pL->signalNorm = guess;
        pVal = pValue( tMu(pL) );
      //  std::cout << "up " << guess << " " << pVal << std::endl;
    }
   
    pL->signalNorm = guess;
    muUp = findtMu90(pL);

    while(guess > 1e-6 &&  pVal > .1)
    {
        guess*=.97;
        pL->signalNorm = guess;
        pVal = pValue( tMu(pL) );
       // std::cout << "down " << guess << " - " << pVal << std::endl;
    }
    if(guess < 1e-6)
    {
        muDown = 0;
    } 
    else
    {    
        pL->signalNorm = guess;      //update guess to find other solution
            
        muDown = muUp;
        //if converge to same solution, try again
        while ( fabs(muDown - muUp)/muUp < 1e-4)
        {
            muDown = findtMu90(pL);
            pL->signalNorm*=0.8;
        }
    }
        
    //print out result
    //std::cout << "90\% confidence interval of mu = " << muDown << " - " << muUp << std::endl << " based on " << pL->detectors[pL->detj].nEvents << " events (" << nEvents << " expected).\n";
    if (muDown > muUp)
    {
        CI[1]=muDown;
        CI[0]=muUp;
    }
    else
    {
        CI[0]=muDown;
        CI[1]=muUp;
    }
    
    return CI;

}

void MCtestConfInterval(paramList *pL)
{

    int detj = 0;
    int total = 0;
    pL->asimov = 0;
    int i;
    double *CI;
    
    for (i=0; i<1000; i++ )
    {
        pL->signalNorm = 1.0;
        pL->detectors[detj].BgNorm = 1.0;
        for(int fluxj=0; fluxj < pL->source.numFlux; fluxj++)
            pL->source.nuFluxNorm[fluxj] = 1.0;
        
        CI = confIntervalSM(pL);    
        if ( CI[0] < 1 && CI[1] > 1 )
            total+=1;
        std::cout << (double)total/(i+1) * 100 << "%\n";
    }
    std::cout << "proportion of MC results within confidence limit: " << (double)total/(i+1) * 100 << "%\n";
}

void confIntVsExposure(paramList *pL)
{
    int detj = 0;
    double *CI;
    char filename[100];
    std::ofstream outfile;
    
    if(pL->elecScat)
        sprintf(filename, "%sCLe_%s_%s_SM.dat",pL->root,pL->detectors[detj].name,pL->source.name);
    else
        sprintf(filename, "%sCLn_%s_%s_SM.dat",pL->root,pL->detectors[detj].name,pL->source.name);
        
    std::cout << "writing output to: " << filename << std::endl;    
    outfile.open(filename,std::ios::out);
    outfile   << "exposure   lower   upper    sigma\n";
        
    double exp = pL->detectors[detj].exposure;
    pL->detectors[detj].exposure = 10 + 90*pL->source.isSolar[0];
    double increment = pow( exp/10, .05);
    double sigma;
    
    while (pL->detectors[detj].exposure <= exp*increment)
    {
        pL->signalNorm = 1.0;
        pL->detectors[detj].BgNorm = 1.0;
        for(int fluxj=0; fluxj < pL->source.numFlux; fluxj++)
            pL->source.nuFluxNorm[fluxj] = 1.0;
        
        CI = confIntervalSM(pL);    
        sigma = sqrt(q0SM(pL));
        
        //print out result
        std::cout << pL->detectors[detj].exposure << "  " << CI[0] << " - " << CI[1] << ", sigma = " << sigma << std::endl;
        outfile   << pL->detectors[detj].exposure << "  " << CI[0] << "  " << CI[1] << "  " << sigma  << std::endl;
        
        pL->detectors[detj].exposure *= increment;
    }
    outfile.close();

    
}


