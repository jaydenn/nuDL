//definition of detector struct
#ifndef STDIO_H
    #include <stdio.h>
#endif
#ifndef GSL_SPLINE_H
    #include <gsl/gsl_spline.h>
#endif

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
    
    gsl_spline *background;
    gsl_interp_accel *accel;
    
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
        nIso=-1; AM=-1; isoA[0]=-1; isoFrac[0]=-1; ErL=0; ErU=-1; bg=-1; eff=-1; res=-1, nEvents=0;
                
        background = gsl_spline_alloc(gsl_interp_linear,500);
        accel = gsl_interp_accel_alloc();
    }
};

