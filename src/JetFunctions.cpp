#include "JetFunctions.h"
#include "integrator.h"

const int frame=1;
double gs, gssq, alphas_med, NPscale, alpha_s_fix, mD2_prefactor, rhoG_prefactor, Tmin;
const double fm2invGeV = 5.068;
const double Nf = 3;
const double Nc = 3;
const double CA = Nc;
const double CF = (Nc*Nc-1)/2/Nc;
const double TwoCFminusCA = 2*CF-CA;
const double dA = 8;
const double dF = 3;
const double Tf = 1./2.;
const double t0 = 0.6*5.068;
const double T0 = 0.5;
const double T3_0 = std::pow(T0, 3);
const double beta0 = 11/3*Nc - 2*Nf/3;
const double FourPioverbeta0 = 4*M_PI/beta0;
const double LambdaQCD_sq_LO = std::pow(0.15, 2);
const double Q2min = 0.25;
const double t_thermal = 1.2*fm2invGeV;

void set_parameters(double _gs, double _screen_scale, double _fix_alphas) {
    gs = _gs;
    gssq = std::pow(gs, 2);
    alphas_med = gssq/4./M_PI;
    NPscale = _screen_scale;
    alpha_s_fix = _fix_alphas;
    mD2_prefactor = (1.+Nf/6.)*gssq;
    rhoG_prefactor = 1.202*3/2/M_PI*alphas_med * (2*dA+CF/CA*3*dF*Nf);
    Tmin = 0.17;
    std::cout << "RhoG prefactor = " << rhoG_prefactor << std::endl;
}

bool GetQns(double x, double kx, double ky, double qx, double qy, double px, double py, 
            double toverx1mxEplus, double mD2,
            std::vector<std::vector<double> > & An, std::vector<double> & Phin, double & Phi13, double & Phi24,
            int type){
    
    double Qn[8][2], Qnsq[8];
    An.clear(); Phin.clear();
    An.resize(8); Phin.resize(8);
    for (auto & it : An) it.resize(2);
    
    double mx = 1.-x;
    double basic_x = x*kx - mx*px, basic_y = x*ky - mx*py;
    
    if(type==0){
        Qn[0][0] = basic_x + mx*qx;        Qn[0][1] = basic_y + mx*qy;
        Qn[1][0] = basic_x;                Qn[1][1] = basic_y;
        Qn[2][0] = basic_x - x*qx;         Qn[2][1] = basic_y - x*qy;
        Qn[3][0] = basic_x - qx;           Qn[3][1] = basic_y - qy;
        for(int i=0;i<4;i++) Qnsq[i] = std::pow(Qn[i][0],2) + std::pow(Qn[i][1],2) + mD2;
        for(int i=0;i<4;i++){
            Phin[i] = 1. - std::cos(Qnsq[i]*toverx1mxEplus);
            for (int j=0;j<2;j++) An[i][j] = Qn[i][j]/Qnsq[i];
        }
        Phi13 = 1. - std::cos((Qnsq[0]-Qnsq[2])*toverx1mxEplus);
        Phi24 = 1. - std::cos((Qnsq[1]-Qnsq[3])*toverx1mxEplus);
    }
    else{
        Qn[4][0] = basic_x + mx*(qx+kx);        Qn[4][1] = basic_y + mx*(qy+ky);
        Qn[5][0] = basic_x + mx*kx;             Qn[5][1] = basic_y + mx*ky;
        Qn[6][0] = basic_x - x*qx + mx*kx;      Qn[6][1] = basic_y - x*qy + mx*ky;
        Qn[7][0] = basic_x + qx + mx*kx;        Qn[7][1] = basic_y + qy + mx*ky;
        for(int i=4;i<8;i++) Qnsq[i] = std::pow(Qn[i][0],2) + std::pow(Qn[i][1],2) + mD2;
        for(int i=4;i<8;i++){
            Phin[i] = 1. - std::cos(Qnsq[i]*toverx1mxEplus);
            for (int j=0;j<2;j++) An[i][j] = Qn[i][j]/Qnsq[i];
        }
    }
    return true;
}

double q2qgJet_WFoverlap_unmeasured(double x, double kx, double ky, double qx, double qy, 
                                  double Eplus, double R, double t, double T){
    double ConeSq = std::pow(std::tan(R/2),2);
    double x1mxEplus = x*(1-x)*Eplus;
    double toverx1mxEplus = t/x1mxEplus;
    double mD2 = GetmD2(T);
    double res = 0.;
    double px = 0., py = 0.;
    double mx = 1. - x;
    std::vector<std::vector<double> > An;
    std::vector<double> Phin;
    double Phi13, Phi24;
    bool status;
    // RR contribution
    if (frame==0) { px = qx - kx; py = qy - ky; }
    else { px = -kx; py = -ky; }
    if (ConeSq < (std::pow(x*kx-mx*px, 2)+std::pow(x*ky-mx*py, 2))/std::pow(x1mxEplus, 2)) {
        status = GetQns(x, kx, ky, qx, qy, px, py, toverx1mxEplus, mD2, An, Phin, Phi13, Phi24, 0);
        double M11 = prodsum(An[1], An[1]), M02 = prodsum(An[0], An[2]),
               M22 = prodsum(An[2], An[2]), M00 = prodsum(An[0], An[0]),
               M12 = prodsum(An[1], An[2]), M01 = prodsum(An[0], An[1]);
        res   += CA*M02*Phi13
               + TwoCFminusCA*( M00 - M01 )*Phin[0]
               + CA*( M22 - M12 )*Phin[2]
               + CA*( M22 - M02 )*Phin[2]
               + CA*( M00 - M02 )*Phin[0];
    }
    
    // RV contribution
    px = - kx; py = - ky; 
    if (ConeSq < (std::pow(x*kx-mx*px, 2)+std::pow(x*ky-mx*py, 2))/std::pow(x1mxEplus, 2)) {
        if (frame == 0) status = GetQns(x, kx, ky, qx, qy, px, py, toverx1mxEplus, mD2, An, Phin, Phi13, Phi24, 0);
        double M11 = prodsum(An[1], An[1]), M13 = prodsum(An[1], An[3]);
        res += - CA*M13*Phi24
               - CA*( M11 - M13 )*Phin[1];
    }
    
    /*
    // VR contribution
    if (frame==0) { px = qx; py = qy; }
    else { px = 0.; py = 0.; }
    status = GetQns(x, kx, ky, qx, qy, px, py, toverx1mxEplus, mD2, An, Phin, Phi13, Phi24, 1);    
    res +=  -  2*CF*prodsum(An[4], An[4])*Phin[4]
            + TwoCFminusCA*prodsum(An[4], An[5])*Phin[4]
            + CA*prodsum(An[5], An[6])*Phin[6];
    
    // VV contribution
    px = 0.; py = 0.;
    if (frame == 0) status = GetQns(x, kx, ky, qx, qy, px, py, toverx1mxEplus, mD2, An, Phin, Phi13, Phi24, 1);
    res +=  + CA*prodsum(An[5], An[7])*Phin[7]
            - CA*prodsum(An[5], An[5])*Phin[5];
    */
    
    return -res;
}



double q2qgJet_WFoverlap_measures_q(double x, double kx, double ky, double qx, double qy, 
                                double Eplus, double R, double t, double T){
    double ConeSq = std::pow(std::tan(R/2),2);
    double x1mxEplus = x*(1-x)*Eplus;
    double toverx1mxEplus = t/x1mxEplus;
    double mD2 = GetmD2(T);
    double res = 0.;
    double px = 0., py = 0.;
    double mx = 1. - x;
    std::vector<std::vector<double> > An;
    std::vector<double> Phin;
    double Phi13, Phi24;
    bool status;
    // RR contribution
    if (frame==0) { px = qx - kx; py = qy - ky; }
    else { px = -kx; py = -ky; }
    if (ConeSq <= (std::pow(x*kx-mx*px, 2)+std::pow(x*ky-mx*py, 2))/std::pow(x1mxEplus, 2)) {
        status = GetQns(x, kx, ky, qx, qy, px, py, toverx1mxEplus, mD2, An, Phin, Phi13, Phi24, 0);
        double M11 = prodsum(An[1], An[1]), M02 = prodsum(An[0], An[2]),
               M22 = prodsum(An[2], An[2]), M00 = prodsum(An[0], An[0]),
               M12 = prodsum(An[1], An[2]), M01 = prodsum(An[0], An[1]);
        res +=   CF*M11 *0
               + CA*M02*Phi13
               + TwoCFminusCA*( M00 - M01 )*Phin[0]
               + CA*( M22 - M12 )*Phin[2]
               + CA*( M22 - M02 )*Phin[2]
               + CA*( M00 - M02 )*Phin[0];
    }
    
    // RV contribution
    px = - kx; py = - ky; 
    if (ConeSq <= (std::pow(x*kx-mx*px, 2)+std::pow(x*ky-mx*py, 2))/std::pow(x1mxEplus, 2)) {
        if (frame == 0) status = GetQns(x, kx, ky, qx, qy, px, py, toverx1mxEplus, mD2, An, Phin, Phi13, Phi24, 0);
        double M11 = prodsum(An[1], An[1]), M13 = prodsum(An[1], An[3]);
        res += - CF*M11*0
               - CA*M13*Phi24
               - CA*( M11 - M13 )*Phin[1];
    }
            
    return res;
}


double g2ggJet_WFoverlap_unmeasured(double x, double kx, double ky, double qx, double qy, 
                                    double Eplus, double R, double t, double T){
    double ConeSq = std::pow(std::tan(R/2),2);
    double x1mxEplus = x*(1-x)*Eplus;
    double toverx1mxEplus = t/x1mxEplus;
    double mD2 = GetmD2(T);
    double res = 0.;
    double px = 0., py = 0.;
    double mx = 1. - x;
    std::vector<std::vector<double> > An;
    std::vector<double> Phin;
    double Phi13, Phi24;
    bool status;
    // only use frame=1
    if (ConeSq < (std::pow(kx, 2)+std::pow(ky, 2))/std::pow(x1mxEplus, 2)) {
        std::vector<double> A{kx, ky}; double A2 = prodsum(A, A) + mD2; A[0]/=A2; A[1]/=A2;
        std::vector<double> B{kx+x*qx, ky+x*qy}; double B2 = prodsum(B, B) + mD2; B[0]/=B2; B[1]/=B2;
        std::vector<double> C{kx-(1-x)*qx, ky-(1-x)*qy}; double C2 = prodsum(C, C) + mD2; C[0]/=C2; C[1]/=C2;
        std::vector<double> D{kx-qx, ky-qy}; double D2 = prodsum(D, D) + mD2; D[0]/=D2; D[1]/=D2;
        double PhiA = 1. - std::cos(A2*toverx1mxEplus);
        double PhiB = 1. - std::cos(B2*toverx1mxEplus);
        double PhiC = 1. - std::cos(C2*toverx1mxEplus);
        double PhiAD = 1. - std::cos((A2-D2)*toverx1mxEplus);
        double PhiBC = 1. - std::cos((B2-C2)*toverx1mxEplus);
        double MAA = prodsum(A, A), MBB = prodsum(B, B), MCC = prodsum(C, C),
               MAB = prodsum(A, B), MAC = prodsum(A, C), MBC = prodsum(B, C), MAD = prodsum(A, D);
        res += CA*( (MBB-MAB)*PhiB + (MBB-MBC)*PhiB + (MCC-MAC)*PhiC )
             + CA*( (MCC-MBC)*PhiC - (MAA-MAD)*PhiA
                   + MBC*PhiBC - MAD*PhiAD ) ;
    }
    return -res;
}

double g2ggJet_WFoverlap_measures_g(double x, double kx, double ky, double qx, double qy, 
                                    double Eplus, double R, double t, double T){
    double ConeSq = std::pow(std::tan(R/2),2);
    double x1mxEplus = x*(1-x)*Eplus;
    double toverx1mxEplus = t/x1mxEplus;
    double mD2 = GetmD2(T);
    double res = 0.;
    double px = 0., py = 0.;
    double mx = 1. - x;
    std::vector<std::vector<double> > An;
    std::vector<double> Phin;
    double Phi13, Phi24;
    bool status;
    // only use frame=1
    if (ConeSq < (std::pow(kx, 2)+std::pow(ky, 2))/std::pow(x1mxEplus, 2)) {
        std::vector<double> A{kx, ky}; double A2 = prodsum(A, A) + mD2; A[0]/=A2; A[1]/=A2;
        std::vector<double> B{kx+x*qx, ky+x*qy}; double B2 = prodsum(B, B) + mD2; B[0]/=B2; B[1]/=B2;
        std::vector<double> C{kx-(1-x)*qx, ky-(1-x)*qy}; double C2 = prodsum(C, C) + mD2; C[0]/=C2; C[1]/=C2;
        std::vector<double> D{kx-qx, ky-qy}; double D2 = prodsum(D, D) + mD2; D[0]/=D2; D[1]/=D2;
        double PhiA = 1. - std::cos(A2*toverx1mxEplus);
        double PhiB = 1. - std::cos(B2*toverx1mxEplus);
        double PhiC = 1. - std::cos(C2*toverx1mxEplus);
        double PhiAD = 1. - std::cos((A2-D2)*toverx1mxEplus);
        double PhiBC = 1. - std::cos((B2-C2)*toverx1mxEplus);
        double MAA = prodsum(A, A), MBB = prodsum(B, B), MCC = prodsum(C, C),
               MAB = prodsum(A, B), MAC = prodsum(A, C), MBC = prodsum(B, C), MAD = prodsum(A, D);
        res += CA*( (MBB-MAB)*PhiB + (MBB-MBC)*PhiB + (MCC-MAC)*PhiC )
             + CA*( (MCC-MBC)*PhiC - (MAA-MAD)*PhiA
                   + MBC*PhiBC - MAD*PhiAD ) ;
    }
            
    return res;
}

double g2qqbarJet_WFoverlap_unmeasured(double x, double kx, double ky, double qx, double qy, 
                                    double Eplus, double R, double t, double T){
    double ConeSq = std::pow(std::tan(R/2),2);
    double x1mxEplus = x*(1-x)*Eplus;
    double toverx1mxEplus = t/x1mxEplus;
    double mD2 = GetmD2(T);
    double res = 0.;
    double px = 0., py = 0.;
    double mx = 1. - x;
    std::vector<std::vector<double> > An;
    std::vector<double> Phin;
    double Phi13, Phi24;
    bool status;
    // only use frame=1
    if (ConeSq < (std::pow(kx, 2)+std::pow(ky, 2))/std::pow(x1mxEplus, 2)) {
        std::vector<double> A{kx, ky}; double A2 = prodsum(A, A) + mD2; A[0]/=A2; A[1]/=A2;
        std::vector<double> B{kx+x*qx, ky+x*qy}; double B2 = prodsum(B, B) + mD2; B[0]/=B2; B[1]/=B2;
        std::vector<double> C{kx-(1-x)*qx, ky-(1-x)*qy}; double C2 = prodsum(C, C) + mD2; C[0]/=C2; C[1]/=C2;
        std::vector<double> D{kx-qx, ky-qy}; double D2 = prodsum(D, D) + mD2; D[0]/=D2; D[1]/=D2;
        double PhiA = 1. - std::cos(A2*toverx1mxEplus);
        double PhiB = 1. - std::cos(B2*toverx1mxEplus);
        double PhiC = 1. - std::cos(C2*toverx1mxEplus);
        double PhiAD = 1. - std::cos((A2-D2)*toverx1mxEplus);
        double PhiBC = 1. - std::cos((B2-C2)*toverx1mxEplus);
        double MAA = prodsum(A, A), MBB = prodsum(B, B), MCC = prodsum(C, C),
               MAB = prodsum(A, B), MAC = prodsum(A, C), MBC = prodsum(B, C), MAD = prodsum(A, D);
        res += CA*(MBB-MAB)*PhiB + TwoCFminusCA*(MBB-MBC)*PhiB + CA*(MCC-MAC)*PhiC
             + TwoCFminusCA*( (MCC-MBC)*PhiC - (MAA-MAD)*PhiA
                             + MBC*PhiBC - MAD*PhiAD ) ;
    }
    return -res;
}

double g2qqbarJet_WFoverlap_measures_q(double x, double kx, double ky, double qx, double qy, 
                                    double Eplus, double R, double t, double T){
    double ConeSq = std::pow(std::tan(R/2),2);
    double x1mxEplus = x*(1-x)*Eplus;
    double toverx1mxEplus = t/x1mxEplus;
    double mD2 = GetmD2(T);
    double res = 0.;
    double px = 0., py = 0.;
    double mx = 1. - x;
    std::vector<std::vector<double> > An;
    std::vector<double> Phin;
    double Phi13, Phi24;
    bool status;
    // only use frame=1
    if (ConeSq < (std::pow(kx, 2)+std::pow(ky, 2))/std::pow(x1mxEplus, 2)) {
        std::vector<double> A{kx, ky}; double A2 = prodsum(A, A) + mD2; A[0]/=A2; A[1]/=A2;
        std::vector<double> B{kx+x*qx, ky+x*qy}; double B2 = prodsum(B, B) + mD2; B[0]/=B2; B[1]/=B2;
        std::vector<double> C{kx-(1-x)*qx, ky-(1-x)*qy}; double C2 = prodsum(C, C) + mD2; C[0]/=C2; C[1]/=C2;
        std::vector<double> D{kx-qx, ky-qy}; double D2 = prodsum(D, D) + mD2; D[0]/=D2; D[1]/=D2;
        double PhiA = 1. - std::cos(A2*toverx1mxEplus);
        double PhiB = 1. - std::cos(B2*toverx1mxEplus);
        double PhiC = 1. - std::cos(C2*toverx1mxEplus);
        double PhiAD = 1. - std::cos((A2-D2)*toverx1mxEplus);
        double PhiBC = 1. - std::cos((B2-C2)*toverx1mxEplus);
        double MAA = prodsum(A, A), MBB = prodsum(B, B), MCC = prodsum(C, C),
               MAB = prodsum(A, B), MAC = prodsum(A, C), MBC = prodsum(B, C), MAD = prodsum(A, D);
        res += CA*(MBB-MAB)*PhiB + TwoCFminusCA*(MBB-MBC)*PhiB + CA*(MCC-MAC)*PhiC
             + TwoCFminusCA*( (MCC-MBC)*PhiC - (MAA-MAD)*PhiA
                             + MBC*PhiBC - MAD*PhiAD ) ;
    }
            
    return res;
}

///////////// With medium average //////////////////////
double Compute_MedAvg_unmeasured_q2qg2J(MediumProfile * Med, TABProfile * TAB, 
                                     double Ejet, double R, double tmax, double xymax, double & err){
    double kmax = Ejet/2;
    double qmax = kmax;
    double Eplus = 2*Ejet;
    double kscale = std::sqrt(Ejet/tmax);
    double lnkmin = std::log(1), lnkmax = std::log(1+kmax/kscale);
    double lnqmin = std::log(1), lnqmax = std::log(1+qmax/kscale);
    double zmin = NPscale/Ejet;
    double zmid = std::max(0.5, 1.0-T3_0*tmax*t0/Ejet);
    double zmax = 1.-NPscale/Ejet;
    double umin = z2u(zmin), umid = z2u(zmid), umax = z2u(zmax);
    double prefactor = 1./2./std::pow(M_PI,2)*alphas_med*CF/M_PI
                       /(2*M_PI * TAB->GetIntegratedTAB());
    
    const int Ndim = 8;
    double xmin1[Ndim]={0.,   umin, lnkmin, lnqmin, 0.,   0., -xymax, -xymax};
    double xmax1[Ndim]={tmax, umid, lnkmax, lnqmax, M_PI, 2*M_PI, xymax, xymax};
    
    double xmin2[Ndim]={0.,   umid, lnkmin, lnqmin, 0.,   0., -xymax, -xymax};
    double xmax2[Ndim]={tmax, umax, lnkmax, lnqmax, M_PI, 2*M_PI, xymax, xymax};
    
    auto dF = [Med, TAB, Ejet, Eplus, R, kscale](double * x_){
        double t = x_[0],  u = x_[1], 
               lk = x_[2], lq = x_[3], 
               phikq = x_[4], 
               phi_jet = x_[5], 
               jet_x0 = x_[6], 
               jet_y0 = x_[7];
        double z = u2z(u); 
        double k = kscale * (std::exp(lk) - 1.);
        double q = kscale * (std::exp(lq) - 1.);
        
        if (k*k > z*(1-z)*Ejet*Ejet/4) return 0.;
        double tab = TAB->GetTAB(jet_y0, jet_x0); // note that the convention of TRNETo is the transpose of that for hydro
        if (tab<1e-5) return 0.0;
        double T = 0.0;
        Med->GetTemp(t, jet_x0+std::cos(phi_jet)*t, jet_y0+std::sin(phi_jet)*t, T);
        if (t<t_thermal) T = T*std::pow(t/t_thermal,1./3.);
        if (T<Tmin) return 0.0;
        double mD2 = GetmD2(T);
        double mD = std::sqrt(mD2);
        if (q*q > Ejet*mD   ) return 0.;
        
        double Jacobian = 2.*M_PI * 2. * k * q * Jacobian_z2u(z) * (k + kscale) * (q + kscale);
        double alpha_s_rad = alpha_s(q*q);
        return tab * alpha_s_rad * rhoG(T) / std::pow(q*q + mD2, 2) * Pqq(z)
             * q2qgJet_WFoverlap_unmeasured(z, k, 0., q*std::cos(phikq), q*std::sin(phikq), Eplus, R, t, T) * Jacobian;
    };
    double err2;
    double res2 = vegas(dF, Ndim, xmin2, xmax2, err2, 200000);
    double err1;
    double res1 = vegas(dF, Ndim, xmin1, xmax1, err1, 200000);
    err = std::sqrt(std::pow(err1,2)+std::pow(err2,2));
    err *= prefactor;
    return (res1+res2)*prefactor;
}


double Compute_MedAvg_measured_q2qg2J(MediumProfile * Med, TABProfile * TAB, 
                                   double z, double Ejet, double R, double tmax, double xymax, double & err){

    double kmax = Ejet/2;
    double qmax = kmax;
    double Eplus = 2*Ejet;
    double kscale = std::sqrt(Ejet/tmax);
    double lnkmin = std::log(1), lnkmax = std::log(1+kmax/kscale);
    double lnqmin = std::log(1), lnqmax = std::log(1+qmax/kscale);
    double prefactor = 1./2./std::pow(M_PI,2)*alphas_med*CF/M_PI
                       /(2*M_PI * TAB->GetIntegratedTAB());
                           
    const int Ndim = 7;
    double xmin[Ndim]={0.,   lnkmin, lnqmin, 0.,   0., -xymax, -xymax};
    double xmax[Ndim]={tmax, lnkmax, lnqmax, M_PI, 2*M_PI, xymax, xymax};
   
        
    auto dF = [Med, TAB, Ejet, Eplus, R, kscale, z](double * x_){
        double t = x_[0], 
               lk = x_[1], lq = x_[2], 
               phikq = x_[3],
               phi_jet = x_[4], 
               jet_x0 = x_[5], 
               jet_y0 = x_[6];
        double k = kscale * (std::exp(lk) - 1.);
        double q = kscale * (std::exp(lq) - 1.);
        
        if (k*k > z*(1-z)*Ejet*Ejet/4) return 0.;
        double tab = TAB->GetTAB(jet_y0, jet_x0); // note that the convention of TRNETo is the transpose of that for hydro
        if (tab<1e-5) return 0.0;
        double T = 0.0;
        Med->GetTemp(t, jet_x0+std::cos(phi_jet)*t, jet_y0+std::sin(phi_jet)*t, T);
        if (t<t_thermal) T = T*std::pow(t/t_thermal,1./3.);
        if (T<Tmin) return 0.0;
        double mD2 = GetmD2(T);
        double mD = std::sqrt(mD2);
        if (q*q > Ejet*mD  
        || z*Ejet<NPscale || (1-z)*Ejet<NPscale 
         ) return 0.;
        
        double Jacobian = 2.*M_PI * 2. * k * q * (k + kscale) * (q + kscale);
        double alpha_s_rad = alpha_s(q*q);
        return tab * alpha_s_rad * rhoG(T) / std::pow(q*q + mD2, 2)
             * q2qgJet_WFoverlap_measures_q(z, k, 0., q*std::cos(phikq), q*std::sin(phikq), Eplus, R, t, T) * Jacobian;
    };
    double res = vegas(dF, Ndim, xmin, xmax, err);
    err *= prefactor*Pqq(z);
    return res*prefactor*Pqq(z);
}



///////////// g2gg
double Compute_MedAvg_unmeasured_g2gg2J(MediumProfile * Med, TABProfile * TAB, 
                                     double Ejet, double R, double tmax, double xymax, double & err){
    double kmax = Ejet/2;
    double qmax = kmax;
    double Eplus = 2*Ejet;
    double kscale = std::sqrt(Ejet/tmax);
    double lnkmin = std::log(1), lnkmax = std::log(1+kmax/kscale);
    double lnqmin = std::log(1), lnqmax = std::log(1+qmax/kscale);
    double zmin = NPscale/Ejet;
    double zmid = std::max(0.5, 1.0-T3_0*tmax*t0/Ejet);
    double zmax = 1.-NPscale/Ejet;
    double umin = z2u(zmin), umid = z2u(zmid), umax = z2u(zmax);
    double prefactor = 1./2./std::pow(M_PI,2)*alphas_med*CF/M_PI
                       /(2*M_PI * TAB->GetIntegratedTAB());
    
    const int Ndim = 8;
    double xmin1[Ndim]={0.,   umin, lnkmin, lnqmin, 0.,   0., -xymax, -xymax};
    double xmax1[Ndim]={tmax, umid, lnkmax, lnqmax, M_PI, 2*M_PI, xymax, xymax};
    
    double xmin2[Ndim]={0.,   umid, lnkmin, lnqmin, 0.,   0., -xymax, -xymax};
    double xmax2[Ndim]={tmax, umax, lnkmax, lnqmax, M_PI, 2*M_PI, xymax, xymax};
    
    auto dF = [Med, TAB, Ejet, Eplus, R, kscale](double * x_){
        double t = x_[0],  u = x_[1], 
               lk = x_[2], lq = x_[3], 
               phikq = x_[4], 
               phi_jet = x_[5], 
               jet_x0 = x_[6], 
               jet_y0 = x_[7];
        double z = u2z(u); 
        double k = kscale * (std::exp(lk) - 1.);
        double q = kscale * (std::exp(lq) - 1.);
        
        if (k*k > z*(1-z)*Ejet*Ejet/4) return 0.;
        double tab = TAB->GetTAB(jet_y0, jet_x0); // note that the convention of TRNETo is the transpose of that for hydro
        if (tab<1e-5) return 0.0;
        double T = 0.0;
        Med->GetTemp(t, jet_x0+std::cos(phi_jet)*t, jet_y0+std::sin(phi_jet)*t, T);
        if (t<t_thermal) T = T*std::pow(t/t_thermal,1./3.);
        if (T<Tmin) return 0.0;
        double mD2 = GetmD2(T);
        double mD = std::sqrt(mD2);
        if (q*q > Ejet*mD  ) return 0.;
        
        double Jacobian = 2.*M_PI * 2. * k * q * Jacobian_z2u(z) * (k + kscale) * (q + kscale);
        double alpha_s_rad = alpha_s(q*q);
        return tab * alpha_s_rad * rhoG(T) / std::pow(q*q + mD2, 2) * Pgg(z)
             * g2ggJet_WFoverlap_unmeasured(z, k, 0., q*std::cos(phikq), q*std::sin(phikq), Eplus, R, t, T) * Jacobian;
    };
    double err2;
    double res2 = vegas(dF, Ndim, xmin2, xmax2, err2, 200000);
    double err1;
    double res1 = vegas(dF, Ndim, xmin1, xmax1, err1, 200000);
    err = std::sqrt(std::pow(err1,2)+std::pow(err2,2));
    err *= prefactor;
    return (res1+res2)*prefactor;
}


double Compute_MedAvg_measured_g2gg2J(MediumProfile * Med, TABProfile * TAB, 
                                   double z, double Ejet, double R, double tmax, double xymax, double & err){

    double kmax = Ejet/2;
    double qmax = kmax;
    double Eplus = 2*Ejet;
    double kscale = std::sqrt(Ejet/tmax);
    double lnkmin = std::log(1), lnkmax = std::log(1+kmax/kscale);
    double lnqmin = std::log(1), lnqmax = std::log(1+qmax/kscale);
    double prefactor = 1./2./std::pow(M_PI,2)*alphas_med*CF/M_PI
                       /(2*M_PI * TAB->GetIntegratedTAB());
                           
    const int Ndim = 7;
    double xmin[Ndim]={0.,   lnkmin, lnqmin, 0.,   0., -xymax, -xymax};
    double xmax[Ndim]={tmax, lnkmax, lnqmax, M_PI, 2*M_PI, xymax, xymax};
   
        
    auto dF = [Med, TAB, Ejet, Eplus, R, kscale, z](double * x_){
        double t = x_[0], 
               lk = x_[1], lq = x_[2], 
               phikq = x_[3],
               phi_jet = x_[4], 
               jet_x0 = x_[5], 
               jet_y0 = x_[6];
        double k = kscale * (std::exp(lk) - 1.);
        double q = kscale * (std::exp(lq) - 1.);
        
        if (k*k > z*(1-z)*Ejet*Ejet/4) return 0.;
        double tab = TAB->GetTAB(jet_y0, jet_x0); // note that the convention of TRNETo is the transpose of that for hydro
        if (tab<1e-5) return 0.0;
        double T = 0.0;
        Med->GetTemp(t, jet_x0+std::cos(phi_jet)*t, jet_y0+std::sin(phi_jet)*t, T);
        if (t<t_thermal) T = T*std::pow(t/t_thermal,1./3.);
        if (T<Tmin) return 0.0;
        double mD2 = GetmD2(T);
        double mD = std::sqrt(mD2);
        if (q*q > Ejet*mD  || z*Ejet<NPscale || (1-z)*Ejet<NPscale ) return 0.;
        
        double Jacobian = 2.*M_PI * 2. * k * q * (k + kscale) * (q + kscale);
        double alpha_s_rad = alpha_s(q*q);
        return tab * alpha_s_rad * rhoG(T) / std::pow(q*q + mD2, 2)
             * g2ggJet_WFoverlap_measures_g(z, k, 0., q*std::cos(phikq), q*std::sin(phikq), Eplus, R, t, T) * Jacobian;
    };
    double res = vegas(dF, Ndim, xmin, xmax, err);
    err *= prefactor*Pgg(z);
    return res*prefactor*Pgg(z);
}



///////////// g2qqbar
double Compute_MedAvg_unmeasured_g2qqbar2J(MediumProfile * Med, TABProfile * TAB, 
                                     double Ejet, double R, double tmax, double xymax, double & err){
    double kmax = Ejet/2;
    double qmax = kmax;
    double Eplus = 2*Ejet;
    double kscale = std::sqrt(Ejet/tmax);
    double lnkmin = std::log(1), lnkmax = std::log(1+kmax/kscale);
    double lnqmin = std::log(1), lnqmax = std::log(1+qmax/kscale);
    double zmin = NPscale/Ejet;
    double zmid = std::max(0.5, 1.0-T3_0*tmax*t0/Ejet);
    double zmax = 1.-NPscale/Ejet;
    double umin = z2u(zmin), umid = z2u(zmid), umax = z2u(zmax);
    double prefactor = 1/2./std::pow(M_PI,2)*alphas_med*CF/M_PI
                       /(2*M_PI * TAB->GetIntegratedTAB());
    
    const int Ndim = 8;
    double xmin1[Ndim]={0.,   umin, lnkmin, lnqmin, 0.,   0., -xymax, -xymax};
    double xmax1[Ndim]={tmax, umid, lnkmax, lnqmax, M_PI, 2*M_PI, xymax, xymax};
    
    double xmin2[Ndim]={0.,   umid, lnkmin, lnqmin, 0.,   0., -xymax, -xymax};
    double xmax2[Ndim]={tmax, umax, lnkmax, lnqmax, M_PI, 2*M_PI, xymax, xymax};
    
    auto dF = [Med, TAB, Ejet, Eplus, R, kscale](double * x_){
        double t = x_[0],  u = x_[1], 
               lk = x_[2], lq = x_[3], 
               phikq = x_[4], 
               phi_jet = x_[5], 
               jet_x0 = x_[6], 
               jet_y0 = x_[7];
        double z = u2z(u); 
        double k = kscale * (std::exp(lk) - 1.);
        double q = kscale * (std::exp(lq) - 1.);
        
        if (k*k > z*(1-z)*Ejet*Ejet/4) return 0.;
        double tab = TAB->GetTAB(jet_y0, jet_x0); // note that the convention of TRNETo is the transpose of that for hydro
        if (tab<1e-5) return 0.0;
        double T = 0.0;
        Med->GetTemp(t, jet_x0+std::cos(phi_jet)*t, jet_y0+std::sin(phi_jet)*t, T);
        if (t<t_thermal) T = T*std::pow(t/t_thermal,1./3.);
        if (T<Tmin) return 0.0;
        double mD2 = GetmD2(T);
        double mD = std::sqrt(mD2);
        if (q*q > Ejet*mD  ) return 0.;
        
        double Jacobian = 2.*M_PI * 2. * k * q * Jacobian_z2u(z) * (k + kscale) * (q + kscale);
        
        double alpha_s_rad = alpha_s(q*q);
        return tab *alpha_s_rad * rhoG(T) / std::pow(q*q + mD2, 2) * Pqg(z)
             * g2qqbarJet_WFoverlap_unmeasured(z, k, 0., q*std::cos(phikq), q*std::sin(phikq), Eplus, R, t, T) * Jacobian;
    };
    double err2;
    double res2 = vegas(dF, Ndim, xmin2, xmax2, err2, 200000);
    double err1;
    double res1 = vegas(dF, Ndim, xmin1, xmax1, err1, 200000);
    err = std::sqrt(std::pow(err1,2)+std::pow(err2,2));
    err *= prefactor;
    return (res1+res2)*prefactor;
}


double Compute_MedAvg_measured_g2qqbar2J(MediumProfile * Med, TABProfile * TAB, 
                                   double z, double Ejet, double R, double tmax, double xymax, double & err){

    double kmax = Ejet/2;
    double qmax = kmax;
    double Eplus = 2*Ejet;
    double kscale = std::sqrt(Ejet/tmax);
    double lnkmin = std::log(1), lnkmax = std::log(1+kmax/kscale);
    double lnqmin = std::log(1), lnqmax = std::log(1+qmax/kscale);
    double prefactor = 1./2./std::pow(M_PI,2)*alphas_med*CF/M_PI
                       /(2*M_PI * TAB->GetIntegratedTAB());
                           
    const int Ndim = 7;
    double xmin[Ndim]={0.,   lnkmin, lnqmin, 0.,   0., -xymax, -xymax};
    double xmax[Ndim]={tmax, lnkmax, lnqmax, M_PI, 2*M_PI, xymax, xymax};
   
        
    auto dF = [Med, TAB, Ejet, Eplus, R, kscale, z](double * x_){
        double t = x_[0], 
               lk = x_[1], lq = x_[2], 
               phikq = x_[3],
               phi_jet = x_[4], 
               jet_x0 = x_[5], 
               jet_y0 = x_[6];
        double k = kscale * (std::exp(lk) - 1.);
        double q = kscale * (std::exp(lq) - 1.);
        
        if (k*k > z*(1-z)*Ejet*Ejet/4) return 0.;
        double tab = TAB->GetTAB(jet_y0, jet_x0); // note that the convention of TRNETo is the transpose of that for hydro
        if (tab<1e-5) return 0.0;
        double T = 0.0;
        Med->GetTemp(t, jet_x0+std::cos(phi_jet)*t, jet_y0+std::sin(phi_jet)*t, T);
        if (t<t_thermal) T = T*std::pow(t/t_thermal,1./3.);
        if (T<Tmin) return 0.0;
        double mD2 = GetmD2(T);
        double mD = std::sqrt(mD2);
        if (q*q > Ejet*mD  || z*Ejet<NPscale || (1-z)*Ejet<NPscale ) return 0.;
        
        double Jacobian = 2.*M_PI * 2. * k * q * (k + kscale) * (q + kscale);
        
        double alpha_s_rad = alpha_s(q*q);
        return tab * alpha_s_rad * rhoG(T) / std::pow(q*q + mD2, 2)
             * g2qqbarJet_WFoverlap_measures_q(z, k, 0., q*std::cos(phikq), q*std::sin(phikq), Eplus, R, t, T) * Jacobian;
    };
    double res = vegas(dF, Ndim, xmin, xmax, err);
    err *= prefactor*Pqg(z);
    return res*prefactor*Pqg(z);
}
















/////////////// Compute collisional energy loss //////////////
// we consider the collisional energy loss is always out of cone, those get back into the jet from medium response is R^2 suppressed for small R jet.
double Compute_MedAvg_Eloss_coll_q(MediumProfile * Med, TABProfile * TAB, 
                                 double Ejet, double tmax, double xymax, double & err){

    const int Ndim = 4;
    double xmin[Ndim]={0.,   0., -xymax, -xymax};
    double xmax[Ndim]={tmax, 2*M_PI, xymax, xymax};
    
    auto dF = [Med, TAB, Ejet](double * x_){
        double t = x_[0],  
               phi_jet = x_[1], 
               jet_x0  = x_[2], 
               jet_y0  = x_[3];
        
        double tab = TAB->GetTAB(jet_y0, jet_x0); // note that the convention of TRNETo is the transpose of that for hydro
        if (tab<1e-5) return 0.0;
        double T = 0.0;
        Med->GetTemp(t, jet_x0+std::cos(phi_jet)*t, jet_y0+std::sin(phi_jet)*t, T);
        if (t<t_thermal) T = T*std::pow(t/t_thermal,1./3.);
        double mD2 = GetmD2(T);
        double mD = std::sqrt(mD2);
        double Q2jet = Ejet*mD;
        if (T<Tmin || Q2jet<mD2) return 0.0;
        double alphas_jet = alpha_s(Q2jet);
        return tab * CF * M_PI * alphas_med * alphas_jet * (1.+Nf/6.) * T*T * std::log(Q2jet/mD2);
    };

    double res = vegas(dF, Ndim, xmin, xmax, err, 200000);
    err /= (2*M_PI * TAB->GetIntegratedTAB());
    return res/(2*M_PI * TAB->GetIntegratedTAB());
}
