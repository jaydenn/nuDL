#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "SMrate.h"
#include "BSMrate.h"
#include "detectorFunctions.h"
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

int calcRates(paramList *pList)
{
    double ErkeV;
    std::string filename = pList->root;    
    std::ofstream outfile;
    std::string BSMname;
    
    switch(pList->BSM)
    {
        case 0:
        {
            BSMname ="SM";
        }
        case 1:
        {
            BSMname = "scalar";
            break;
        }
        case 2:
        {
            BSMname = "pseudoS";
            break;
        }
        case 3:
        {
            BSMname = "vector";
            break;
        }
        case 4:
        {
            BSMname = "axial";
            break;
        }
        case 5:
        {
            BSMname = "sterile";
            break;
        }
    }
    
    //format output streams
    std::cout << std::setiosflags(std::ios::scientific) << std::setprecision(4);
    outfile   << std::setiosflags(std::ios::scientific) << std::setprecision(5);
    
    //output model
    int masskeV = (int)(pList->mMed*1e6);
    if (pList->BSM != 0 && pList->BSM != 5)
        std::cout << "\nBSM rate for " << masskeV << " keV "<< BSMname << " mediator\n";
      
    for(int detj=0; detj < pList->ndet; detj++)
    {
        //open file
        if(pList->nucScat)
            filename += "Rate_NR_";
        else
            filename += "Rate_ER_";
        
        filename += BSMname + pList->detectors[detj].name + pList->source.name;
        outfile.open(filename,std::ios::out);
        
        if( !outfile )
        {
            std::cout << "output file could not be created (does the directory in specified root path exist?" << std::endl;
            return 1;
        }
        
        std::cout << "------------------------\n";
        std::cout << detj+1 << ". " << pList->detectors[detj].name << std::endl;
        std::cout << "------------------------\n";
        std::cout << "  Total rates: \n"; 
        std::cout << "     SM  = " << intSMrate(  pList->detectors[detj].ErL, pList->detectors[detj].ErU, pList, detj)         << " events/kg/day" << std::endl;
        std::cout << "     BG  = " << intBgRate(  pList->detectors[detj], pList->detectors[detj].ErL, pList->detectors[0].ErU) << " events/kg/day" << std::endl;
        if( pList->BSM != 0 )
            std::cout << "     BSM = " << intBSMrate( pList->detectors[detj].ErL, pList->detectors[detj].ErU, pList, detj,1)       << " events/kg/day" << std::endl;
        
        std::cout << "  Differential rates: \n"; 
        std::cout << "    Er (keV)        SM dR/dE        BG dR/dE"; 
        if( pList->BSM != 0 )
            std::cout << "        BSM dN/dE";
        std::cout << "(events/kg/day/keV)\n";
        
        outfile   << "    Er (keV)         SM dR/dE         BG dR/dE";
        if( pList->BSM != 0 )
            outfile << "         BSM dN/dE";
        outfile << "(events/kg/day/keV)\n";
        
        int skip=0;            
        for (int i=0; i<501; i+=1)
        {
            if(pList->logBins == 1)
                ErkeV = pow(10, log10(pList->detectors[detj].ErL) + (double)i*(log10(pList->detectors[detj].ErU)-log10(pList->detectors[detj].ErL))/500);
            else
                ErkeV = pList->detectors[detj].ErL + (double)i*(pList->detectors[detj].ErU-pList->detectors[detj].ErL)/500;

            if( skip++ % 5 == 0 ) //only print out every fifth value to terminal
            {
                std::cout << "    " << ErkeV << "      " << diffSMrate( ErkeV, pList, detj) << "      " << diffBgRate( pList->detectors[detj], ErkeV);
                if (pList->BSM != 0)
                    std::cout << "      " << diffBSMrate( ErkeV, pList, detj, 1);
                std::cout << std::endl;
            }
            
            outfile   << "    " << ErkeV << "      " << diffSMrate( ErkeV, pList, detj) << "      " << diffBgRate( pList->detectors[detj], ErkeV);
            if (pList->BSM != 0)
                outfile << "      " << diffBSMrate( ErkeV, pList, detj, 1);
            outfile << std::endl;
        }
        
        outfile.close();
    }
   
}

int calcRatesThreshold(paramList *pList)
{
    double ErkeV;
    std::string filename = pList->root;    
    std::ofstream outfile;
    std::string BSMname;
    
    switch(pList->BSM)
    {
        case 0:
        {
            BSMname ="SM";
        }
        case 1:
        {
            BSMname = "scalar";
            break;
        }
        case 2:
        {
            BSMname = "pseudoS";
            break;
        }
        case 3:
        {
            BSMname = "vector";
            break;
        }
        case 4:
        {
            BSMname = "axial";
            break;
        }
        case 5:
        {
            BSMname = "sterile";
            break;
        }
    }
    
    //format output streams
    std::cout << std::setiosflags(std::ios::scientific) << std::setprecision(4);
    outfile   << std::setiosflags(std::ios::scientific) << std::setprecision(5);
    
    //output model
    int masskeV = (int)(pList->mMed*1e6);
    if (pList->BSM != 0 && pList->BSM != 5)
        std::cout << "\nBSM rate for " << masskeV << " keV "<< BSMname << " mediator\n";
      
    for(int detj=0; detj < pList->ndet; detj++)
    {
        //open file
        if(pList->nucScat)
            filename += "RateTh_NR_";
        else
            filename += "RateTh_ER_";
        
        filename += BSMname + pList->detectors[detj].name + pList->source.name;
        outfile.open(filename,std::ios::out);
        
        if( !outfile )
        {
            std::cout << "output file could not be created (try creating the directory \'results\'" << std::endl;
            return 1;
        }
        
        std::cout << "------------------------\n";
        std::cout << detj+1 << ". " << pList->detectors[detj].name << std::endl;
        std::cout << "------------------------\n";
        std::cout << "  total rates above threshold: \n"; 
        std::cout << "    ErTh (keV)      SM rate         BG rate";
        if(pList->BSM != 0)
            std::cout << "        BSM rate";
        std::cout <<" (events/kg/day)\n";
        outfile   << "    ErTh (keV)      SM rate         BG rate";
        if(pList->BSM != 0)
            outfile   << "         BSM rate"; 
        outfile <<" (events/kg/day)\n";

        int skip=0;            
        for (int i=0; i<501; i+=1)
        {
            if(pList->logBins == 1)
                ErkeV = pow(10, log10(pList->detectors[detj].ErL) + (double)i*(log10(pList->detectors[detj].ErU)-log10(pList->detectors[detj].ErL))/500);
            else
                ErkeV = pList->detectors[detj].ErL + (double)i*(pList->detectors[detj].ErU-pList->detectors[detj].ErL)/500;

            if( skip++ % 5 == 0 ) //only print out every fifth value to terminal
            {            
                std::cout << "    " << ErkeV << "      " << intSMrate( ErkeV, pList->detectors[detj].ErU, pList, detj) << "      " << intBgRate( pList->detectors[detj], ErkeV, pList->detectors[detj].ErU);
                if (pList->BSM != 0)
                    std::cout << "      " << intBSMrate( ErkeV, pList->detectors[detj].ErU, pList, detj, 1);
                std::cout << std::endl;
            }
            
            outfile   << "    " << ErkeV << "      " << intSMrate( ErkeV, pList->detectors[detj].ErU, pList, detj) << "      " << intBgRate( pList->detectors[detj], ErkeV, pList->detectors[detj].ErU);
             if (pList->BSM != 0)
                outfile << "      " << intBSMrate( ErkeV, pList->detectors[detj].ErU, pList, detj, 1);
            outfile << std::endl; 
        }
        
        outfile.close();
    }
    
}    
