#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>
#include <H5Cpp.h>
#include "JetFunctions.h"
#include "MediumLoader.h"




int main(int argc, char * argv[]){

    // input  
    double sqrts = atof(argv[1]); 
    double R = atof(argv[2]);
    std::string label(argv[3]);
    double gs_med = atof(argv[4]); 
    double NP_scale = atof(argv[5]);
    double alphas_fix = atof(argv[6]);
    std::string hydrofolder = std::string(argv[7]);
    std::string outputfolder = std::string(argv[8]);
    
    set_parameters(gs_med, NP_scale, alphas_fix);
    
    double res, err, err1, err2;
    int NE = 21, Nz=101;
    double Ejet0=8, Ejet1=sqrts/2;
    
    std::stringstream name0, name1, name2;
    name0 << "Coll_eloss-"<<label<<".dat";
    name1 << "Unmeasured-"<<label<<".dat";
    name2 << "Measured-"<<label<<".dat";
    
    // Load hydro file
    MediumProfile * Med = new MediumProfile();
    
    Med->ReadHydroGridFile(hydrofolder+"/JetData.h5");

    // Load TAB file
    TABProfile * TAB = new TABProfile();
    TAB->ReadTABGridFile(hydrofolder+"/ic.h5");


    // maximum of QGP grid
    double tmax = 15.0*fm2invGeV;
    double xymax = 10.0*fm2invGeV;
    
    // collisional energy loss
    std::ofstream f0(outputfolder+name0.str());
    for (int i=0; i<NE; i++){
        double s = i*1./(NE-1);
        double Ejet = std::exp(std::log(Ejet0)+s*std::log(Ejet1/Ejet0));
        res = Compute_MedAvg_Eloss_coll_q(Med, TAB, Ejet, tmax, xymax, err);
        std::cout << Ejet << " " << res << " +/- "<< err<<std::endl;
        f0 << Ejet << " " << res << " "<< err<<std::endl; 
    }    
    
    // unresolved jet function
    std::ofstream f1(outputfolder+name1.str());
    for (int i=0; i<NE; i++){
        double s = i*1./(NE-1);
        double Ejet = std::exp(std::log(Ejet0)+s*std::log(Ejet1/Ejet0));
        double resq2qg = Compute_MedAvg_unmeasured_q2qg2J(Med, TAB, Ejet, R, tmax, xymax, err);
        double resg2gg = Compute_MedAvg_unmeasured_g2gg2J(Med, TAB, Ejet, R, tmax, xymax, err);
        double resg2qqbar = Compute_MedAvg_unmeasured_g2qqbar2J(Med, TAB, Ejet, R, tmax, xymax, err);
        std::cout << Ejet << " " << resq2qg << " " << resg2gg << " " << resg2qqbar << std::endl;
        f1 << Ejet << " " << resq2qg << " " << resg2gg << " " << resg2qqbar << std::endl; 
    }
    
    // resolved jet function
    std::ofstream f2(outputfolder+name2.str());
    for (int i=0; i<NE; i++){
        double s = i*1./(NE-1);
        double Ejet = std::exp(std::log(Ejet0)+s*std::log(Ejet1/Ejet0));
        for (int j=0; j<Nz; j++){
            double u = -9+18*j*1.0/(Nz-1);
            double z = u2z(u);
            double resq2qg = Compute_MedAvg_measured_q2qg2J(Med, TAB, z, Ejet, R, tmax, xymax, err1)
                           + Compute_MedAvg_measured_q2qg2J(Med, TAB, 1-z, Ejet, R, tmax, xymax, err2) ;
            double resg2gg = Compute_MedAvg_measured_g2gg2J(Med, TAB, z, Ejet, R, tmax, xymax, err1)
                           + Compute_MedAvg_measured_g2gg2J(Med, TAB, 1-z, Ejet, R, tmax, xymax, err2) ;
            double resg2qqbar = Compute_MedAvg_measured_g2qqbar2J(Med, TAB, z, Ejet, R, tmax, xymax, err1)
                              + Compute_MedAvg_measured_g2qqbar2J(Med, TAB, 1-z, Ejet, R, tmax, xymax, err2) ;
                    
            std::cout << Ejet 
                      << " "   << std::setprecision(5) << z 
                      << " "   << std::setprecision(5) << z*(1-z)*resq2qg 
                      << " "   << std::setprecision(5) << z*(1-z)*resg2gg
                      << " "   << std::setprecision(5) << z*(1-z)*resg2qqbar
                      << std::endl;
            f2 << Ejet 
               << " "   << std::setprecision(9) << z 
               << " "   << std::setprecision(9) << z*(1-z)*resq2qg 
               << " "   << std::setprecision(9) << z*(1-z)*resg2gg
               << " "   << std::setprecision(9) << z*(1-z)*resg2qqbar
               << std::endl;
        }
    }
    
}
