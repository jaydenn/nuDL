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
#include "nuRate.h"
#include "SMrate.h"
#include "BSMrate.h"


double my_LmuM(const gsl_vector *v, void *params)
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
    
    l -= logLikelihood(pL);
}

double my_Lmu(const gsl_vector *v, void *params)
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
    
    l -= logLikelihood(pL);
    
    return l;

}

double findMaxLmuM(paramList *pL)
{
    size_t iter = 0;
    int status;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_multimin_function my_func;
    
    my_func.f = my_LmuM;
    my_func.params = (void *)pL;

    //start point and step size
    gsl_vector *x,*dx;
    my_func.n = 2 + + pL->source.numFlux;

    x = gsl_vector_alloc ( my_func.n );
    dx = gsl_vector_alloc ( my_func.n );
    gsl_vector_set (x, 0, pL->signalNorm);
    gsl_vector_set(dx, 0, pL->signalNorm/10);
    
    for(int i=1; i < my_func.n; i++)
    {
        gsl_vector_set (x, i, 1.0);
        gsl_vector_set(dx, i, .05);
    }
    
    gsl_multimin_fminimizer_set (s, &my_func, x, dx);

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
    
    return -LS;
}


double findMaxLmu(paramList *pL)
{
    size_t iter = 0;
    int status;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_multimin_function my_func;

    
    my_func.f = my_Lmu;
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
  //      cout << "       " << iter << "  c " << pL->w.coeff << " " <<  gsl_vector_get (s->x, 0) << " " <<  gsl_vector_get (s->x, 1) << " " <<  gsl_vector_get (s->x, 2) << " " << s->fval << " " << gsl_multimin_fminimizer_size(s) << endl; 
    }
    while (iter < 900 && gsl_multimin_fminimizer_size(s)>.0001);///s->fval > 1e-2 && !status);
    if(iter==900)
        std::cout << "non-convergence size = " << gsl_multimin_fminimizer_size(s)/s->fval  << " > .0001  " <<  std::endl;
    
    double Lmu = s->fval;
    
    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free (dx);
            
    return -Lmu;
}

//statistic for discovery
double qMu(paramList *pL)
{	 
    double maxLmu = findMaxLmu( pL );
    double q;

    if( pL->signalNorm < 0)
        q = 0.0;
    else
        q = - 2 * ( maxLmu - pL->maxL );  //no -ve because functions return -loglike

    if (q < 0 )
    {   
        if (q < -1e-5)
            std::cout << "warning maximum in likelihood not reached max L = " << pL->maxL << " < " << maxLmu << std::endl;
        return 0;
    }
    else
        return q;
}

double my_qMu(const gsl_vector *v, void *params)
{

    double simSeed = 0;
    paramList *pL = (paramList *)params;
    pL->signalNorm  = gsl_vector_get(v, 0);
    
    return pow( sqrt(qMu(pL)) - 1.282, 2);  //90% confidence exclusion
    
}


double findqMu90(paramList *pL)
{

    size_t iter = 0;
    int status;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;

    gsl_vector *x,*dx;
    gsl_multimin_function my_func;

    my_func.n = 1;
    my_func.f = my_qMu;
    my_func.params = (void *)pL;

    //only need to calc max likelihood once per mass
    //pL->maxL = findMaxLmuM( pL );
 
    //cout << "maxL " << pL->p.maxL << endl;
    //cout << "maxL " << pL->w.coeff << endl;
    
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

void exclusionLimit(paramList *pL, int detj)
{

    char filename[100];
    std::ofstream outfile;
    
    if(pL->elecScat)
        sprintf(filename, "%sexclusionE_%c%c_BSM%d.dat",pL->root,pL->detectors[0].name[0],pL->detectors[0].name[1],pL->BSM);
    else
        sprintf(filename, "%sexclusionN_%c%c_BSM%d.dat",pL->root,pL->detectors[0].name[0],pL->detectors[0].name[1],pL->BSM);
    
    outfile.open(filename,std::ios::out);
    std::cout << "writing output to: " << filename << std::endl;
    
    //find a good starting point
    double mu  = 1e-4;
    double BSM = intBSMrate( pL->detectors[detj].ErL, pL->detectors[detj].ErU, pL, detj, mu);
    double BG  = intBgRate(pL->detectors[detj], pL->detectors[detj].ErL, pL->detectors[detj].ErU);         
    double SM  = intSMrate( pL->detectors[detj].ErL, pL->detectors[detj].ErU, pL, detj);
    
    starting_Ex_guess_loop:
        while(fabs(BSM) < SM/10 )
        {
            mu*=1.02;
            BSM = intBSMrate( pL->detectors[detj].ErL, pL->detectors[detj].ErU, pL, detj, mu);
           // std::cout << SM << " " << BG << " " << BSM << " " << mu << std::endl;
        }
    
        //produce a null result
        pL->signalNorm = 0;
        generateBinnedData( pL, pL->detj, 1, 0);
        pL->maxL = logLikelihood(pL);
        double q=10;
        while(mu > 1e-8 && sqrt(q) > 1.28)
        {
            mu*=.9;
            pL->signalNorm = mu;    
            q = qMu( pL );
            //std::cout << mu << " " << q << std::endl;
        }
  
//    std::cout << mu << "  " << pL->maxL << std::endl;      
    
    double coup;
    pL->signalNorm = mu;
    pL->detj = detj;
    
    while (pL->mMed < 1)
    {

        mu = findqMu90(pL);
        
        if (mu==mu) //check for NAN
        {
            if(pL->BSM==3 || pL->BSM==4)
                coup=mu*pL->C;
            else
                coup=sqrt(mu)*pL->C;
                
            //print out result
            std::cout << pL->mMed << "  " << coup << std::endl;
            outfile   << pL->mMed << "  " << coup << std::endl;
            
            pL->signalNorm = mu;      //update guess
        }
        else
        {
            mu = 1e-6;
            goto starting_Ex_guess_loop;
        }
        
        for(int i=0; i < pL->source.numFlux; i++)
            pL->source.nuFluxNorm[i] = 1;
        pL->mMed*=1.2; //increment mass
        
        //reinitialize BSM rates
        for(int fluxj=0; fluxj< pL->source.numFlux; fluxj++)
        {
            pL->SMinterference1=1;  pL->SMinterference2=0;
            rateInit( pL, detj, fluxj, &BSMrate,  pL->detectors[detj].signalBSM1[fluxj]);
	        if(pL->BSM==3 || pL->BSM==4)
	        {
	                pL->SMinterference2=1; pL->SMinterference1=0;
	                rateInit( pL, detj, fluxj, &BSMrate,  pL->detectors[detj].signalBSM2[fluxj]);
	        }
	        pL->SMinterference1=pL->SMinterference2=1;
        }    
    }
    outfile.close();

}

