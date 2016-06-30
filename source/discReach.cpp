#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_multimin.h>
#include "likelihood.h"
#include "monteCarlo.h"
#include "detectors.h"
#include "assert.h"
#include <string.h>

double bestFitOP(parameterList *pL, int OP);

double my_LS(const gsl_vector *v, void *params)
{

    parameterList *pL = (parameterList *)params;
    double l;

    for(int i=0; i<pL->ndet; i++)
    {
        pL->detectors[i].normB8  = gsl_vector_get(v, 3);
    }
    pL->w.coeff = gsl_vector_get(v, 4);
    
    l = - logLikelihood( pL, 1) + 19.5 *( pow(pL->detectors[0].normB8 - 1,2) ) ;    
    return l;
}

double my_L0(const gsl_vector *v, void *params)
{

    parameterList *pL = (parameterList *)params;
    pL->w.coeff = 0;
    double l;

    for(int i=0; i<pL->ndet; i++)
    {
        pL->detectors[i].normB8  = gsl_vector_get(v, 3);
    }
       
    l = - logLikelihood( pL, 1) + 19.5 * pow(pL->detectors[0].normB8 - 1,2) ;
    return l;
}

double findMaxLS(parameterList *pL)
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
    my_func.n = 5;

    x = gsl_vector_alloc (5);
    dx = gsl_vector_alloc (5);
    gsl_vector_set (x, 0, 1.0);
    gsl_vector_set(dx, 0, .05);
    
    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, 5);

    gsl_multimin_fminimizer_set (s, &my_func, x, dx);


    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate (s);       
        //cout << "       " <<iter << " " <<  gsl_vector_get (s->x, 0) << " " <<  gsl_vector_get (s->x, 1) << " " <<  gsl_vector_get (s->x, 2) << " " <<  gsl_vector_get (s->x, 3) << " " << s->fval << endl; 
    }
    while (iter < 2000 && gsl_multimin_fminimizer_size(s)>.001); //s->fval>1e-2);
    
   if(iter==2000)
        cout << "LS non-convergence size = " << gsl_multimin_fminimizer_size(s) << " > .001  " << endl;
    
    double LS =  s->fval;

    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free (dx);
    
    return LS;
}


double findMaxL0(parameterList *pL)
{

    size_t iter = 0;
    int status;

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_multimin_function my_func;

    
    my_func.f = my_L0;
    my_func.params = (void *)pL;

    //start point
    gsl_vector *x,*dx;
    my_func.n = 4;
    
    x = gsl_vector_alloc (4);
    dx = gsl_vector_alloc (4);
    gsl_vector_set (x, 0, 1.0); 
    gsl_vector_set(dx, 0, .05);

    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, 4);

    gsl_multimin_fminimizer_set (s, &my_func, x, dx);

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate (s);
        //cout << "       " << iter << " " <<  gsl_vector_get (s->x, 0) << " " <<  gsl_vector_get (s->x, 1) << " " <<  gsl_vector_get (s->x, 2) << " " << s->fval << endl; 
    }
    while (iter < 2000 && gsl_multimin_fminimizer_size(s)>.001);///s->fval > 1e-2 && !status);
    if(iter==2000)
        cout << "L0 non-convergence size = " << gsl_multimin_fminimizer_size(s)  << " > .001  " <<  endl;
    
    double L0 = s->fval;
    
    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free (dx);
            
    return L0;
}

//statistic for discovery
double q0(parameterList *pL)
{	 

    double maxL0 = findMaxL0( pL );
    double maxL  = findMaxLS( pL );
    
    return 2 * ( maxL0 - maxL );  //no -ve because functions return -loglike
 
}

double my_q0(const gsl_vector *v, void *params)
{

    double simSeed = 0;
    parameterList *pL = (parameterList *)params;    
   
    pL->w.coeff = gsl_vector_get(v, 0);       
    pL->p.coeffn[pL->OP][0] =  gsl_vector_get(v, 0);
    generateBinnedData( &(pL->w), &(pL->detectors[i]), 1, simSeed);
 
    return pow( sqrt(q0(pL)) - 4.28, 2);  //arbitrary function with a minima at 3 sigma
}


double findCoeff3sig(parameterList *pL)
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
    
    gsl_vector_set (x, 0, pL->p.coeffn[pL->OP][0]);
    gsl_vector_set(dx, 0, pL->p.coeffn[pL->OP][1]);

    T = gsl_multimin_fminimizer_nmsimplex2;
    s = gsl_multimin_fminimizer_alloc (T, 1);

    gsl_multimin_fminimizer_set (s, &my_func, x, dx);

    double prevC =0;
    double convCount=0;
    do
    {

        prevC = gsl_vector_get(s->x,0);
        status = gsl_multimin_fminimizer_iterate (s);
        if(prevC==gsl_vector_get(s->x,0))
        {
            convCount++;
            if(convCount>100)
            {
                convCount=0;
                gsl_vector_set(x,0,1.1*gsl_vector_get(s->x,0));
                gsl_multimin_fminimizer_free (s);
                s = gsl_multimin_fminimizer_alloc (T, 1);   
                gsl_multimin_fminimizer_set (s, &my_func, x, dx);
            }
            
        }
        iter++;
        //cout << iter << " " << gsl_vector_get(s->x,0) << ", q' = " << s->fval << ", size " << gsl_multimin_fminimizer_size(s) << endl;
       
    }
    while (iter < 2000 && s->fval > .0002 && !status); //under 1% error in 4.28 sigma value
        
    double c = fabs(gsl_vector_get(s->x, 0));
    
    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free (dx);

    if(iter==2000)
    {
        cout << "non-convergence f = " << s->fval << " > .0002" << endl;
        return -1;
    }
    else
        return c;
}

double find_C (const gsl_vector *v, void *params)
{
    parameterList *pL = (parameterList *)params;

    pL->w.coeff = gsl_vector_get (v, 0);
    if (gsl_vector_get(v, 0) < 0)
        pL->w.sign = -1;
    else
        pL->w.sign = 1;
        
    return -logLikelihood( pL, 0); 
}

void discLimit(parameterList *pL)
{

    char filename[100];
    std::ofstream outfile;
    
    sprintf(filename, "%snuFloor_%c%c_C%d_%c_low.dat",pL->root,pL->detectors[0].name[0],pL->detectors[0].name[1],pL->OP,N);
    
    cout << "writing output to: " << filename << endl;    
    outfile.open(filename,ios::out);
    
    double mu;
    //InitializeNullHyp();
    //InitializeBSM();

    while (//some scanning param )
    {

        coeff = findCoeff3sig(pL);

        //print out result
        cout    << pL->w.Mx << "  " << coeff << endl;
        outfile << pL->w.Mx << "  " << coeff << endl;
        
        pL->p.coeffn[pL->OP][0] = mu;      //update guess
        pL->p.coeffn[pL->OP][1] = mu/5;   //update stepsize

        pL->w.Mx *= 2;                 //increment threshold?
            
        InitializeWIMP(pL);
        
    }
    outfile.close();

}

/*
void discEvolution(parameterList *pL)
{
    
    char filename[100];
    std::ofstream outfile;
    char N;
    switch(pL->p.isv)
    {
        case 1:
            N='N';
            break;
        case 2:
            N='p';
            break;
        case 3:
            N='n';
            break;
    }
    if(pL->detectors[0].highLow==0)
        sprintf(filename,"%sdiscEvo_%c%c_C%d_%c_low.dat",pL->root,pL->detectors[0].name[0],pL->detectors[0].name[1],pL->OP,N);
    else
        sprintf(filename,"%sdiscEvo_%c%c_C%d_%c_high.dat",pL->root,pL->detectors[0].name[0],pL->detectors[0].name[1],pL->OP,N);

    outfile.open(filename,ios::out);
    
    if( strcmp(pL->detectors[0].name,"XENON_low") == 0 )
    {
        double mOP[] = {0,6,0,4.22,6,4.83,4.2,6.3,6.3,4.6,4.6,4.6,4.6,4.2,4.8,4.08};
        pL->w.Mx = mOP[pL->OP];
    }
    if( strcmp(pL->detectors[0].name,"XENON_high") == 0 )
    {
        if(pL->p.isv==2)
        {
            double mOP[] = {0,150,0,36,54,48,38,54,130,38,45,48,42,36,39,35};
            pL->w.Mx = mOP[pL->OP];           
        }
        else
	    {
	       double mOP[] = {0,2700,0,38,71,44,35,1e6,87,62,39,51,46,37,67,35};
           pL->w.Mx = mOP[pL->OP];
        }
     //   pL->w.Mx = 100;              
    }
    if( strcmp(pL->detectors[0].name,"GERMANIUM_low") == 0 )
    {
        double mOP[] = {0,6,0,4.22,6,4.83,4.2,6.3,6.3,4.6,4.6,4.6,4.6,4.2,4.8,4.08};
        pL->w.Mx = mOP[pL->OP];
    }
    if( strcmp(pL->detectors[0].name,"GERMANIUM_high") == 0 )
    {
        double mOP[] = {0,1000,0,48,150,70,47,183,131,90,64,78,68,46,100,41};
        pL->w.Mx = mOP[pL->OP];
    }
    if( strcmp(pL->detectors[0].name,"SILICON_low") == 0 )
    {
        double mOP[] = {0,7.6,0,5.2,7.7,5.8,5.3,8.2,8.0,5.6,5.6,5.5,5.5,5.8,5.8,5.2};
        pL->w.Mx = mOP[pL->OP];
    }
    if( strcmp(pL->detectors[0].name,"SILICON_high") == 0 )   
        pL->w.Mx = 50;
    if( strcmp(pL->detectors[0].name,"FLOURINE_low") == 0 )
    {
        double mOP[] = {0,9.1,0,7.0,9.1,7.4,7.0,9.7,9.4,7.4,7.4,7.2,7.4,7.2,7.4,6.6};
        pL->w.Mx = mOP[pL->OP];
    }
    if( strcmp(pL->detectors[0].name,"FLOURINE_high") == 0 )
        pL->w.Mx = 50;
    
    pL->detectors[0].exposure = 1*pL->detectors[0].highLow + 1e-6;
    
    cout << "starting O" << pL->OP << endl;
    InitializeWIMP(pL);
    
    //find a good starting point
    double coeff = bestFitC(pL); 
    pL->p.coeffn[pL->OP][0] = coeff;
    pL->p.coeffn[pL->OP][1] = coeff/5;
    
    double nuEvents;
    cout    << setprecision(4);
    outfile << setprecision(6);
    while (pL->detectors[0].exposure < (1e3 + pL->detectors[0].highLow * 1e8))
    {
        
        coeff = findCoeff3sig(pL);
        if (coeff > 0)
        {
            nuEvents =  pL->detectors[0].exposure * intBgRate(&(pL->detectors[0]),pL->detectors[0].ErL,pL->detectors[0].ErU);
            
            //print out results
            cout    <<  pL->detectors[0].exposure << " " << nuEvents << " " << fabs(coeff) << endl;
            outfile <<  pL->detectors[0].exposure << " " << nuEvents << " " << fabs(coeff) << endl;
            
            pL->p.coeffn[pL->OP][0] = coeff;      //update guess
            pL->p.coeffn[pL->OP][1] = coeff/20;   //update stepsize
        }

        pL->detectors[0].exposure *= 1.05;
                    
    }
    outfile.close();
} */

