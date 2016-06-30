#include <cmath>
#include <cstdio>

double detEff(double Er, int type)
{
    switch( type ) 
    {
        case 0: 
            return 1;
        case 1: 
            return .3;
        case 2: 
            return .4;
        case 3: 
            return .5;
        case 4: 
            return 1;
        default:
            printf("invalid detector efficiency\n"); 
            return NAN;
    }
}

double detBackground(double Er, int type)
{
    switch( type ) 
    {
        case 0: 
            return 1e-99;
        case 1: 
            return 8e-6*(1-.9998);
        case 2: 
            return 4e-6;
        case 3: 
            return 1e-5;
        case 4: 
            return 1e-6;
        default:
            printf("invalid detector background\n"); 
            return NAN; 
    }
}
double detRes(double Er, int type)
{
    switch( type ) 
    {
        case 0: 
            return 0;
        case 1: 
            return 0;
        case 2: 
            return 0;
        case 3: 
            return 0;
        case 4: 
            return 0;
        default:
            printf("invalid detector resolution\n"); 
            return NAN; 
    }
}

