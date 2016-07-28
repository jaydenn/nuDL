#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_multimin.h>
#include <string.h>
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
#include "SMrate.h"
#include "BSMrate.h"


double my_LS(const gsl_vector *v, void *params)
{

    paramList *pL = (paramList *)params;
    double l;

    pL->nuFluxNorm = fabs(gsl_vector_get(v, 0));
    pL->signalNorm = gsl_vector_get(v, 1);
    pL->detectors[pL->detj].BgNorm = fabs(gsl_vector_get(v, 2));
        
    l = - logLikelihood(pL) 
        + 0.5/pow(pL->nuFluxUn,2) * pow(pL->nuFluxNorm - 1,2) 
        + 0.5/pow(pL->detectors[pL->detj].BgUn,2) * pow(pL->detectors[pL->detj].BgNorm - 1,2) ;   
    return l;

}

double my_L0(const gsl_vector *v, void *params)
{

    paramList *pL = (paramList *)params;
    double l;

    pL->nuFluxNorm = gsl_vector_get(v, 0);
    pL->detectors[pL->detj].BgNorm  = gsl_vector_get(v, 1);
        
    l = - logLikelihood(pL) 
        + 0.5/pow(pL->nuFluxUn,2) * pow(pL->nuFluxNorm - 1,2) 
        + 0.5/pow(pL->detectors[pL->detj].BgUn,2) * pow(pL->detectors[pL->detj].BgNorm - 1,2) ;
    return l;

}

double findMaxLS(paramList *pL)
{

    size_t iter = 0;
    int status;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_multimin_function my_func;
    
    my_func.f = my_LS;
    my_func.params = (void *)pL;

    //start point and step size
    gsl_vector *x,*dx;
    my_func.n = 3;

    x = gsl_vector_alloc (3);
    dx = gsl_vector_alloc (3);
    gsl_vector_set (x, 0, 1.0);
    gsl_vector_set(dx, 0, .05);
    gsl_vector_set (x, 1, pL->signalNorm);
    gsl_vector_set(dx, 1, pL->signalNorm/10);
    gsl_vector_set (x, 2, 1.0);
    gsl_vector_set(dx, 2, .05);
    
    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, 3);

    gsl_multimin_fminimizer_set (s, &my_func, x, dx);

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate (s);       
        //std::cout << "       " <<iter << " " <<  gsl_vector_get (s->x, 0) << " " <<  gsl_vector_get (s->x, 1) << " " << s->fval << std::endl; 
    }
    while (iter < 2000 && gsl_multimin_fminimizer_size(s)>.0005); //s->fval>1e-2);
    
   if(iter==2000)
        std::cout << "LS non-convergence size = " << gsl_multimin_fminimizer_size(s) << " > .0005  " << std::endl;
    
    double LS =  s->fval;
    pL->signalNorm = gsl_vector_get(s->x, 1);
    
    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free (dx);
    
    return LS;
}


double findMaxL0(paramList *pL)
{

    size_t iter = 0;
    int status;
    
    double mu = pL->signalNorm; //save current mu for later
    pL->signalNorm = 0;
    
    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_multimin_function my_func;

    my_func.f = my_L0;
    my_func.params = (void *)pL;

    //start point
    gsl_vector *x,*dx;
    my_func.n = 2;
    
    x = gsl_vector_alloc (2);
    dx = gsl_vector_alloc (2);
    gsl_vector_set (x, 0, 1.0);
    gsl_vector_set(dx, 0, .05);
    gsl_vector_set (x, 1, 1.0); 
    gsl_vector_set(dx, 1, .05);

    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, 2);

    gsl_multimin_fminimizer_set (s, &my_func, x, dx);

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate (s);
        //std::cout << "       " << iter << " " <<  gsl_vector_get (s->x, 0) << " " << gsl_vector_get (s->x, 1) << " " << s->fval << std::endl; 
    }
    while (iter < 2000 && gsl_multimin_fminimizer_size(s)>.0005);  //s->fval > 1e-2 && !status);
    if(iter==2000)
        std::cout << "L0 non-convergence size = " << gsl_multimin_fminimizer_size(s)  << " > .0005  " <<  std::endl;
    
    double L0 = s->fval;
    
    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free (dx);
    
    pL->signalNorm = mu;
    
    return L0;
}

//statistic for discovery
double q0(paramList *pL)
{	 

    double maxL0 = -findMaxL0( pL );
    double maxL  = -findMaxLS( pL ); //-ve because functions return -loglike for minimization
    
    if( pL->signalNorm > 0 )
        return - 2 * ( maxL0 - maxL );  
    else
        return 0;
}

double my_q0(const gsl_vector *v, void *params)
{

    double simSeed = 0;
    paramList *pL = (paramList *)params;    
    
    pL->nuFluxNorm = 1;
    pL->detectors[pL->detj].BgNorm = 1;
    pL->signalNorm = gsl_vector_get(v, 0);
    
    generateBinnedData( pL, pL->detj, 1, simSeed);

    return pow( sqrt(q0(pL)) - 4.28, 2);  //arbitrary function with a minima at 3 sigma 90% of the time
}


double findCoeff3sig(paramList *pL)
{

    size_t iter = 0;
    int status;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;

    gsl_vector *x,*dx;
    gsl_multimin_function my_func;

    my_func.n = 1;
    my_func.f = my_q0;
    my_func.params = (void *)pL;

    //start point
    x = gsl_vector_alloc (1);
    dx = gsl_vector_alloc (1);

    gsl_vector_set (x, 0, pL->signalNorm);
    gsl_vector_set(dx, 0, pL->signalNorm/5);

    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, 1);

    gsl_multimin_fminimizer_set (s, &my_func, x, dx);

   pL->printPars();
    do
    {
        status = gsl_multimin_fminimizer_iterate (s);
        iter++;
       // std::cout <<  pL->detectors[0].exposure << "  " << iter << " " << gsl_vector_get(s->x,0) << ", q' = " << s->fval << ", size " << gsl_multimin_fminimizer_size(s) << std::endl;
    }
    while (iter < 2000 && s->fval > .0015 && !status); //under 1% error in 4.28 sigma value
        
    double mu = gsl_vector_get(s->x, 0);
    
    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free (dx);

    if(iter==2000)
    {
        std::cout << "non-convergence f = " << s->fval << " > .0015" << std::endl;
        return NAN;
    }
    else
        return mu;
}

void discLimitEvolution(paramList *pL, int detj)
{

    std::cout << "Starting disc. evolution calculations..." << std::endl;

    char filename[100];
    std::ofstream outfile;
    
    if(pL->elecScat)
        sprintf(filename, "%sdiscEvoE_%c%c_BSM%d.dat",pL->root,pL->detectors[0].name[0],pL->detectors[0].name[1],pL->BSM);
    else
        sprintf(filename, "%sdiscEvoN_%c%c_BSM%d.dat",pL->root,pL->detectors[0].name[0],pL->detectors[0].name[1],pL->BSM);
        
    std::cout << "writing output to: " << filename << std::endl;    
    outfile.open(filename,std::ios::out);
    
         
    //determine first guess for mu, need a mu that gives BSM ~ SM/100
    double BSM = intBSMrate( pL->detectors[detj].ErL, pL->detectors[detj].ErU, pL, detj, 1);
    double BG  = intBgRate(pL->detectors[detj], pL->detectors[detj].ErL, pL->detectors[detj].ErU);         
    double SM  = intSMrate( pL->detectors[detj].ErL, pL->detectors[detj].ErU, pL, detj);
    double mu  = (SM+BG)/(30*BSM);
    
    pL->detj = detj;
    pL->signalNorm = mu;         
    
    while (pL->detectors[detj].exposure < 1e6)
    {

        mu = findCoeff3sig(pL);

        if (mu==mu) //check for NAN
        {
            //print out result
            std::cout << pL->detectors[detj].exposure << "  " << mu*pL->C << std::endl;
            outfile   << pL->detectors[detj].exposure << "  " << mu*pL->C << std::endl;
            
            pL->signalNorm = mu;      //update guess
        }
        
        pL->detectors[detj].exposure*=1.2; //increment exposure
        
    }
    outfile.close();

}

void discLimitVsMmed(paramList *pL, int detj)
{

    std::cout << "Starting disc. limit calculations..." << std::endl;

    char filename[100];
    std::ofstream outfile;
    
    if(pL->elecScat)
        sprintf(filename, "%sDLe_%c%c_BSM%d.dat",pL->root,pL->detectors[0].name[0],pL->detectors[0].name[1],pL->BSM);
    else
        sprintf(filename, "%sDLn_%c%c_BSM%d.dat",pL->root,pL->detectors[0].name[0],pL->detectors[0].name[1],pL->BSM);
        
    std::cout << "writing output to: " << filename << std::endl;    
    outfile.open(filename,std::ios::out);
    
    pL->mMed = 1e-5;
    
    //determine first guess for mu, need a mu that gives BSM ~ SM/100
    double BSM = intBSMrate( pL->detectors[detj].ErL, pL->detectors[detj].ErU, pL, detj, 1);
    double BG  = intBgRate(pL->detectors[detj], pL->detectors[detj].ErL, pL->detectors[detj].ErU);         
    double SM  = intSMrate( pL->detectors[detj].ErL, pL->detectors[detj].ErU, pL, detj);
    double mu  = (SM+BG)/(30*BSM);
        
    pL->signalNorm = mu;
    pL->detj = detj;
    
    while (pL->mMed < 1e-1)
    {

        mu = findCoeff3sig(pL);

        if (mu==mu) //check for NAN
        {
            //print out result
            std::cout << pL->mMed << "  " << mu*pL->C << std::endl;
            outfile   << pL->mMed << "  " << mu*pL->C << std::endl;
            
            pL->signalNorm = mu;      //update guess
        }
        
        pL->mMed*=1.2; //increment exposure
        
    }
    outfile.close();

}


