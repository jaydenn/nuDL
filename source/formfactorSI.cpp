#include <cmath>

#include "physicalConstants.h"

//nuclear form factor as a function of recoil energy
double ffactorSI(double A, double Er) 
{
    //Standard Helm       
    return 9.05322e6 * exp(-1.9394e-5 * A * Er ) * (-0.00692 * sqrt(2.17706 + pow(-0.6 + 1.23 * pow(A,0.3333),2)) * sqrt(A * Er) * cos(0.00692 * sqrt(2.17706 + pow(-0.6 + 1.23 * pow(A,0.3333),2)) * sqrt(A * Er)) + sin( 0.00692 * sqrt(2.17706 + pow(-0.6 + 1.23 * pow(A,0.3333),2)) * sqrt(A * Er) )) / ( pow( 2.17706 + pow(-0.6 + 1.23 * pow(A,0.3333),2),1.5) * pow(A * Er,1.5));
}

double ffactorSIhelm(double A, double q)
{
    double s = 0.9; //in fm
    double r = sqrt( pow(1.23 * pow(A,1./3.) - 0.6,2) + 7./3.*pow(M_PI,2)*pow(0.52,2)-5*pow(s,2)); //in fm
    double qr = q*r/HBARCgevfm;
    
    return 3*( sin(qr) - qr*cos(qr))/pow(qr,3) * exp(-pow(q*s/HBARCgevfm,2)/2.);
}
