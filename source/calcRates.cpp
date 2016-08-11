#include <iostream>
#include <fstream>
#include <iomanip>
#include "SMrate.h"
#include "BSMrate.h"
#include "detectorFunctions.h"
#ifndef PARAMETERSTRUCT_H
    #include "parameterStruct.h"
#endif

int calcRates(paramList *pList)
{
    double ErkeV;
    char filename[60];    
    std::ofstream outfile;
    char BSMname[15];
    
    switch(pList->BSM)
    {
        case 1:
        {
            sprintf(BSMname,"scalar");
            break;
        }
        case 2:
        {
            sprintf(BSMname,"pseudoS");
            break;
        }
        case 3:
        {
            sprintf(BSMname,"vector");
            break;
        }
        case 4:
        {
            sprintf(BSMname,"axial");
            break;
        }
    }
    
    //format output streams
    std::cout << std::setiosflags(std::ios::scientific) << std::setprecision(4);
    outfile   << std::setiosflags(std::ios::scientific) << std::setprecision(5);
    
    //output model
    int masskeV = (int)(pList->mMed*1e6);
    std::cout << "\nBSM rate for " << masskeV << " keV "<< BSMname << " mediator\n";
      
    for(int detj=0; detj < pList->ndet; detj++)
    {
        //open file
        if(pList->nucScat)
            sprintf(filename,"%s%sRateN_%s_%c%c_m%dkeV.dat", pList->root, BSMname, pList->detectors[detj].name, pList->source.name[0], pList->source.name[1],masskeV);
        else
            sprintf(filename,"%s%sRateE_%s_%c%c_m%dkeV.dat", pList->root, BSMname, pList->detectors[detj].name, pList->source.name[0], pList->source.name[1], masskeV);
        outfile.open(filename,std::ios::out);
        
        if(outfile==NULL)
        {
            std::cout << "output file could not be created (try creating the directory \'results\'" << std::endl;
            return 1;
        }
        
        std::cout << "------------------------\n";
        std::cout << detj+1 << ". " << pList->detectors[detj].name << std::endl;
        std::cout << "------------------------\n";
        std::cout << "  total rates: \n"; 
        std::cout << "     SM  = " << intSMrate(  pList->detectors[detj].ErL, pList->detectors[detj].ErU, pList, detj)        << " events/kg/day" << std::endl;
        std::cout << "     BG  = " << intBgRate(  pList->detectors[detj], pList->detectors[detj].ErL, pList->detectors[0].ErU) << " events/kg/day" << std::endl;
        std::cout << "     BSM = " << intBSMrate( pList->detectors[detj].ErL, pList->detectors[detj].ErU, pList, detj,1)        << " events/kg/day" << std::endl;
        
        std::cout << "  differential rates: \n"; 
        std::cout << "    Er (keV)        " << "SM dN/dE        " << "BG dN/dE";
        outfile   << "    Er (keV)         " << "SM dN/dE         " << "BG dN/dE";
        if(pList->BSM!=0)
        {
            std::cout << "        BSM dN/dE (events/kg/day/keV)" << std::endl;
            outfile   << "         BSM dN/dE (events/kg/day/keV)" << std::endl;
        }
        else
        {
            std::cout << " (events/kg/day/keV)" << std::endl;
            outfile   << " (events/kg/day/keV)" << std::endl;
        }
            
        for (int i=0; i<201; i+=1)
        {
            ErkeV = pList->detectors[detj].ErL + (double)i*(pList->detectors[detj].ErU-pList->detectors[detj].ErL)/200;
            std::cout << "    " << ErkeV << "      " << diffSMrate( ErkeV, pList, detj) << "      " << diffBgRate( pList->detectors[detj], ErkeV);
            outfile   << "    " << ErkeV << "      " << diffSMrate( ErkeV, pList, detj) << "      " << diffBgRate( pList->detectors[detj], ErkeV);
            
            if(pList->BSM!=0)
            {
                std::cout << "      " << diffBSMrate( ErkeV, pList, detj, 1) << std::endl; 
                outfile   << "      " << diffBSMrate( ErkeV, pList, detj, 1) << std::endl; 
            }
            else
            {
                std::cout << std::endl;
                outfile   << std::endl;
            }
        }
        
        outfile.close();
    }
    
}    
