#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string> 
#include <stdlib.h>
#include <sstream>
#include "detectorFunctions.h"
#include "nuFlux.h"
#include "SMrate.h"
#include "BSMrate.h"
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif
#ifndef DETECTORSTRUCT_H
	#include "detectorStruct.h"
#endif	

//gets sampling parameters from file
int interactiveInput(paramList *pL) 
{

    std::string stdinstr;
    std::string::size_type sz;
     
    int mode;
    char *ret;
    char temp[400];

    //Mode switch
    std::cout << "Operating modes:\n";
    std::cout << " 1: print out expected rates\n";    
    std::cout << " 2: calculate discovery limits\n";
    choose_a_valid_mode_loop:
        std::cout << "Enter mode: ";
        std::getline(std::cin,stdinstr);
        //std::cin >> mode;     
        mode = atoi(stdinstr.c_str());
        
        if(!(mode==2 || mode==1) || std::cin.fail())
        {
            std::cout << "invalid mode\n";
            goto choose_a_valid_mode_loop;
        }
        
    // root for output files
    std::string root;
    choose_a_valid_file_loop:
        std::cout << "\nEnter root for output filenames (default is results/CN_): ";
        std::getline(std::cin,root);

        if(std::cin.fail())
        {
            std::cout << "invalid mode\n";
            goto choose_a_valid_file_loop;
        }
        if(root=="")
        {
            std::cout << "using default \"results/CN_\"\n";
            sprintf(pL->root, "results/CN_");		 
        }
        else
            strcpy(pL->root, root.c_str());		 

    //nuclear or electron scattering?
    char scat;
    std::cout << "\nDo electron (E) or nuclear (N) scattering? ";
    choose_a_valid_scattering_loop:
        std::getline(std::cin,stdinstr);
        scat = stdinstr.c_str()[0];
        
        if(!(scat=='E' || scat=='N'))
        {
            std::cout << "\nplease enter E or N: ";
            goto choose_a_valid_scattering_loop;
        }
        else if(scat=='E')
        {
            pL->elecScat=1;
            pL->nucScat=0;
        }
        else if(scat == 'N' )
        {
            pL->elecScat=0;
            pL->nucScat=1;
        }		 		 

    if( (pL->nucScat == 0 && pL->elecScat == 0) || (pL->nucScat == 1 && pL->elecScat == 1))
    {
        std::cout << "Must choose one of electron or nuclear scattering" << std::endl;
        return -1;
    }
    
    //which BSM model to consider
    std::cout << "\nChoose the BSM interaction to consider:\n";
    std::cout << " 1: scalar\n";    
    std::cout << " 2: pseudoscalar\n";
    std::cout << " 3: vector\n";
    std::cout << " 4: axialvector\n";
    choose_a_valid_BSM_loop:
        std::cout << "Enter 1-4: ";
        std::getline(std::cin,stdinstr);
        //std::cin >> mode;     
        pL->BSM = atoi(stdinstr.c_str());
        
        if(!(pL->BSM==1 || pL->BSM==2 || pL->BSM==3 || pL->BSM==4) || std::cin.fail())
        {
            std::cout << "invalid BSM\n";
            goto choose_a_valid_BSM_loop;
        }
        
    //mediator mass
    choose_a_valid_mMed_loop:
        std::cout << "\nEnter mediator mass in GeV: ";
        std::getline(std::cin,stdinstr);
        pL->mMed = atof(stdinstr.c_str());		 
        
        if(std::cin.fail() || pL->mMed < 0 || stdinstr=="")
        {
            std::cout << "invalid mass\n";
            goto choose_a_valid_mMed_loop;
        }
    
    //source setup
    std::cout << "\nAvailable sources:\n";
    std::ifstream sourceFile("sources.ini");
    
    int i=1;
    std::string sources[10];
    while(!sourceFile.eof())
    {
        std::getline(sourceFile, stdinstr); 
        if(stdinstr[0]=='-' && i < 10)
        {
            std::getline(sourceFile, stdinstr);
            sources[i] = stdinstr;
            std::cout << i << ". " << sources[i] << std::endl;
            i++;
        }
        if(i==10)
            std::cout << "max number of sources reached (10)\n";
    }
    
    choose_a_valid_source_loop:
        std::cout << "\nchoose source from above list: ";
        std::getline(std::cin,stdinstr);
        int s = atof(stdinstr.c_str());		 
        if(std::cin.fail() || s < 1)
        {
            std::cout << "invalid source\n";
            goto choose_a_valid_source_loop;
        }
        
    //initialize reactor flux
    int err = nuFluxInit(pL, sources[s]);
    if (err == 1)
    {
        std::cout << "problem with flux initialization\n";
        return -1;
    }
    
    //Detector setup
    char name[20];
    int ndet;
    std::cout << "\nHow many detectors would you like to use?\n";
    choose_a_valid_ndet_loop:
        std::cout << "Enter num. det. <= 10: ";
        std::getline(std::cin,stdinstr);
        ndet = atoi(stdinstr.c_str());
        
        if( ndet>10 || ndet<1 || std::cin.fail() || stdinstr=="")
        {
            std::cout << "invalid number of detectors\n";
            goto choose_a_valid_ndet_loop;
        }
    
    double exp;
    std::cout << "\nAvailable detectors:\n";
    FILE *detsINI;
    detsINI = fopen("detectors.ini","r");
    if(detsINI==NULL)
    {
	    printf("unable to open detectors.ini\n");
	    return 1;
    }

    i=1;
    char dets[20][100];
    while(!feof(detsINI))
    {
        ret = fgets(temp,200,detsINI);
        if(temp[0]=='#')
        {
            sscanf(temp,"%*c%s",dets[i]);
            std::cout << i << ". " << dets[i] << std::endl;
            i++;
        }
    }
	    
    int j=1;
    int det;
    double dist;
    while(ndet>0)
    {
        std::cout << "\nChoose detector " << j << " from list: ";
        choose_a_valid_det_loop:
            std::getline(std::cin,stdinstr);
            det = atoi(stdinstr.c_str());
            
            if( det>i || det<1 || std::cin.fail() || stdinstr=="")
            {
                std::cout << "\nchoose a valid detector: ";
                goto choose_a_valid_det_loop;
            }
            else
                sprintf(name,"#%s",dets[det]);
        
        if(mode!=1)
        {
            std::cout << "Enter exposure (in kg.days) for " << name << " detector: ";
            choose_a_valid_exp_loop:
                std::getline(std::cin,stdinstr);
                exp = atof(stdinstr.c_str());
                
                if( exp < 0 || std::cin.fail())
                {
                    std::cout << "\nenter a valid exposure: ";
                    goto choose_a_valid_exp_loop;
                }
        }
        else
        {
            exp=1;
        }   
         std::cout << "Enter distance (in m) for " << name << " detector: ";
         choose_a_valid_dist_loop:
            std::getline(std::cin,stdinstr);
            dist = atof(stdinstr.c_str());
            
            if( dist < 0 || std::cin.fail())
            {
                std::cout << "\nenter a valid distance: ";
                goto choose_a_valid_dist_loop;
            }
	    if(newDetector(pL, name, exp, dist)) { std::cout << "Could not create detector, exiting" << std::endl; return -1; }

        ndet--; j++;
    }    
        
    //asimov or random sim?
    std::cout << "\nChoose statistical mode:\n";
    std::cout << " 1: Monte-carlo (random data)\n";    
    std::cout << " 2: Asimov (expected/median data)\n" << std::endl;
    choose_a_valid_MC_loop:
        std::cout << "Enter 1 or 2: ";
        std::getline(std::cin,stdinstr);
        pL->asimov = atoi(stdinstr.c_str())-1;
        
        if(!(pL->asimov==0 || pL->asimov==1) || std::cin.fail() || stdinstr=="")
        {
            std::cout << "invalid mode\n";
            goto choose_a_valid_MC_loop;
        }
    std::cout << std::endl;
    
    return mode;
}

