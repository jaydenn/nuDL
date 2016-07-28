//definition of detector struct
#ifndef STDIO_H
    #include <stdio.h>
#endif
#ifndef GSL_SPLINE_H
    #include <gsl/gsl_spline.h>
#endif
#define DETECTORSTRUCT_H

struct detector {
    char name[20];
    double exposure;
    double AM;
    int nIso;
    int isoZ[10];
    int isoA[10];
    double isoFrac[10];
    double isoSZ[10];
    double isoSN[10];
    double ErL;
    double ErU;
    
    int bg;
    int eff;
    int res;
    
    double *binnedData;     //array of i bins with binnedData[i] number of events per bin
    int nbins;
    double binW;
    double *unbinnedData;   //array of i events which occured at energy unbinnedData[i]
    double nEvents;
    
    double BgNorm, BgUn;  //norm factor for background and fractional uncertainty
    
    gsl_spline *background;
    gsl_interp_accel *accelBg;
    gsl_spline *signalSM;
    gsl_interp_accel *accelSM;
    gsl_spline *signalBSM1;
    gsl_interp_accel *accelBSM1;
    gsl_spline *signalBSM2;
    gsl_interp_accel *accelBSM2;
   
    
    void printDetSpecs()
    {
        printf(" %s\n",name);
        
        for (int i=0; i<nIso; i++)
        {
            printf("    Isotope %d - %2.1f%%\n",i+1,isoFrac[i]*100);
            printf("     A  = %d\n",isoA[i]);    
        }
        
        printf("   %.1f < Er < %.1f\n",ErL,ErU);
        printf("   bg  = %d\n",bg);
        printf("   eff = %d\n",eff);
        printf("   res = %d\n",res);
    }
    
    detector()
    {
        nIso=0; AM=-1; isoA[0]=-1; isoFrac[0]=-1; ErL=0; ErU=-1; bg=-1; eff=-1; res=-1, nEvents=0;
        BgNorm = 1; BgUn = 1e-99;
                
        background = gsl_spline_alloc(gsl_interp_linear,1000);
        accelBg = gsl_interp_accel_alloc();
        signalSM = gsl_spline_alloc(gsl_interp_linear,1000);
        accelSM = gsl_interp_accel_alloc();
        signalBSM1 = gsl_spline_alloc(gsl_interp_linear,1000);
        accelBSM1 = gsl_interp_accel_alloc();
        signalBSM2 = gsl_spline_alloc(gsl_interp_linear,1000);
        accelBSM2 = gsl_interp_accel_alloc();
    
    }
    
    
};

