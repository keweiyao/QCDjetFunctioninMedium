#ifndef JETFUNCTIONS_H
#define JETFUNCTIONS_H

#include <cmath>
#include <vector>
#include "MediumLoader.h"

extern const int frame;
extern double gs, gssq, alphas_med, NPscale, alpha_s_fix, mD2_prefactor, rhoG_prefactor, Tmin; // tunable parameters
extern const double Nf;
extern const double Nc;
extern const double CA;
extern const double CF;
extern const double TwoCFminusCA;
extern const double dA;
extern const double dF;
extern const double Tf;
extern const double t0;
extern const double T0;
extern const double T3_0;
extern const double beta0;
extern const double FourPioverbeta0;
extern const double LambdaQCD_sq_LO;
extern const double Q2min;
extern const double fm2invGeV;
extern const double t_thermal;

void set_parameters(double _gs, double _screen_scale, double _afix);


inline double alpha_s(double Q2){
    if (alpha_s_fix > 0.) return alpha_s_fix;
    else {
        if (Q2<Q2min) return FourPioverbeta0/std::log(Q2min/LambdaQCD_sq_LO);
        return FourPioverbeta0/std::log(Q2/LambdaQCD_sq_LO);
    }
}

inline double GetT3(double t){
    if (t>t0) return T3_0*t0/t;
    else return T3_0*std::pow(t/t0,2);
}

inline double GetT(double t){
    return std::pow(GetT3(t), 1./3.);
}

inline double GetmD2(double T){
    return mD2_prefactor*T*T;
}

inline double rhoG(double T){
    return rhoG_prefactor*std::pow(T,3);
}

inline double prodsum(std::vector<double> & a, std::vector<double> & b){
    return a[0]*b[0] + a[1]*b[1];
}

inline double u2z(double u){
    return 1./(1.+std::exp(-u));
}

inline double z2u(double z){
    return std::log(z/(1.-z));
}

inline double Jacobian_z2u(double z){
    return z*(1.-z);
}

inline double Pqq(double x){
    return (1.+x*x)/(1.-x);
}

inline double Pgq(double x){
    return (1+std::pow(1.-x,2))/x;
}

inline double Pgg(double x){
    return (1+std::pow(1.-x,4)+std::pow(x,4))/(x*(1-x));
}

inline double Pqg(double x){
    return std::pow(x,2)+std::pow(1-x,2);
}


bool GetQns(double x, double kx, double ky, double qx, double qy, double px, double py, 
            double toverx1mxEplus, double mD2,
            std::vector<std::vector<double> > & An, std::vector<double> & Phin, double & Phi13, double & Phi24,
            int type);

// q2qg initiated jets
double q2qgJet_WFoverlap_unmeasured(double x, double kx, double ky, double qx, double qy, 
                                    double Eplus, double R, double t, double T);

double q2qgJet_WFoverlap_measures_q(double x, double kx, double ky, double qx, double qy, 
                                    double Eplus, double R, double t, double T);                        

// g2gg initiated jets
double g2ggJet_WFoverlap_unmeasured(double x, double kx, double ky, double qx, double qy, 
                                    double Eplus, double R, double t, double T);
 
double g2ggJet_WFoverlap_measures_g(double x, double kx, double ky, double qx, double qy, 
                                    double Eplus, double R, double t, double T);

// g2gg initiated jets
double g2qqbarJet_WFoverlap_unmeasured(double x, double kx, double ky, double qx, double qy, 
                                    double Eplus, double R, double t, double T);

double g2qqbarJet_WFoverlap_measures_q(double x, double kx, double ky, double qx, double qy, 
                                    double Eplus, double R, double t, double T);
                                    
// average jet function in medium
double Compute_MedAvg_unmeasured_q2qg2J(MediumProfile * Med, TABProfile* TAB, double Ejet, double R, double tmax, double xymax, double & err);        
double Compute_MedAvg_measured_q2qg2J(MediumProfile * Med, TABProfile* TAB, double z, double Ejet, double R, double tmax, double xymax, double & err);      

double Compute_MedAvg_unmeasured_g2gg2J(MediumProfile * Med, TABProfile* TAB, double Ejet, double R, double tmax, double xymax, double & err);        
double Compute_MedAvg_measured_g2gg2J(MediumProfile * Med, TABProfile* TAB, double z, double Ejet, double R, double tmax, double xymax, double & err);      

double Compute_MedAvg_unmeasured_g2qqbar2J(MediumProfile * Med, TABProfile* TAB, double Ejet, double R, double tmax, double xymax, double & err);        
double Compute_MedAvg_measured_g2qqbar2J(MediumProfile * Med, TABProfile* TAB, double z, double Ejet, double R, double tmax, double xymax, double & err);      


double Compute_MedAvg_Eloss_coll_q(MediumProfile * Med, TABProfile * TAB, double Ejet, double tmax, double xymax, double & err);
                                                     
#endif

