#include "complextimerungekuttafehlberg.h"

//------------------------------------------------------------------------------
ComplexTimeRungeKuttaFehlberg::ComplexTimeRungeKuttaFehlberg(Config *cfg):
    ComplexTimePropagation(cfg)
{
    try{
        epsilon = cfg->lookup("ComplexTimeIntegration.rungeKuttaFehlberg.epsilon");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "ComplexTimeRungeKuttaFehlberg::ComplexTimeRungeKuttaFehlberg(Congif *cfg)"
             << "::Error reading parameter from config object." << endl;
        exit(EXIT_FAILURE);
    }
    i = cx_double(1,0);
}
//------------------------------------------------------------------------------
bool ComplexTimeRungeKuttaFehlberg::stepForward()
{
    cx_vec k1, k2, k3, k4,k5, k6;
    cx_mat m1, m2, m3, m4, m5, m6;
    cx_vec A_, A_1;
    cx_mat C_, C_1;

    bool accepted = false;

    // Computing Runge-Kutta weights
    V->computeNewElements(C);
    h->computeNewElements(C);

    k1 = -i*dt*slater->computeRightHandSideComplexTime(A);
    m1 = -i*dt*orbital->computeRightHandSide(C, A);

    V->computeNewElements(C + m1/4.0);
    h->computeNewElements(C + m1/4.0);

    k2 = -i*dt*slater->computeRightHandSideComplexTime(A + k1/4.0);
    m2 = -i*dt*orbital->computeRightHandSide(C + m1/4.0, A + k1/4.0);

    V->computeNewElements(C + 3.0/32*m1 + 9.0/32*m2);
    h->computeNewElements(C + 3.0/32*m1 + 9.0/32*m2);

    k3 = -i*dt*slater->computeRightHandSideComplexTime(A + 3.0/32*k1 + 9.0/32*k2);
    m3 = -i*dt*orbital->computeRightHandSide(C + 3.0/32*m1 + 9.0/32*m2, A + 3.0/32*k1 + 9.0/32*k2);

    V->computeNewElements(C + 1932.0/2197*m1 - 7200.0/2197*m2 + 7296.0/2197*m3);
    h->computeNewElements(C + 1932.0/2197*m1 - 7200.0/2197*m2 + 7296.0/2197*m3);

    k4 = -i*dt*slater->computeRightHandSideComplexTime(A + 1932.0/2197*k1 - 7200.0/2197*k2 + 7296.0/2197*k3);
    m4 = -i*dt*orbital->computeRightHandSide(C + 1932.0/2197*m1 - 7200.0/2197*m2 + 7296.0/2197*m3,
                                             A + 1932.0/2197*k1 - 7200.0/2197*k2 + 7296.0/2197*k3);

    V->computeNewElements(C + 439.0/216*m1 - 8.0*m2 + 3680.0/513*m3 - 845.0/4104*m4);
    h->computeNewElements(C + 439.0/216*m1 - 8.0*m2 + 3680.0/513*m3 - 845.0/4104*m4);

    k5 = -i*dt*slater->computeRightHandSideComplexTime(A + 439.0/216*k1 - 8.0*k2 + 3680.0/513*k3 - 845.0/4104*k4);
    m5 = -i*dt*orbital->computeRightHandSide(C + 439.0/216*m1 - 8.0*m2 + 3680.0/513*m3 - 845.0/4104*m4,
                                             A + 439.0/216*k1 - 8.0*k2 + 3680.0/513*k3 - 845.0/4104*k4);

    V->computeNewElements(C - 8.0/27*m1 + 2.0*m2 - 3544.0/2565*m3 + 1859.0/4104*m4 - 11.0/40*m5);
    h->computeNewElements(C - 8.0/27*m1 + 2.0*m2 - 3544.0/2565*m3 + 1859.0/4104*m4 - 11.0/40*m5);

    k6 = -i*dt*slater->computeRightHandSideComplexTime(A - 8.0/27*k1 + 2.0*k2 - 3544.0/2565*k3 + 1859.0/4104*k4 - 11.0/40*k5);
    m6 = -i*dt*orbital->computeRightHandSide(C - 8.0/27*m1 + 2.0*m2 - 3544.0/2565*m3 + 1859.0/4104*m4 - 11.0/40*m5,
                                             A - 8.0/27*k1 + 2.0*k2 - 3544.0/2565*k3 + 1859.0/4104*k4 - 11.0/40*k5);


    // Computing new states
    A_  = A + 16.0/135*k1 + 6656.0/12825*k3 + 28561.0/56430*k4 - 9.0/50*k5 + 2.0/55*k6;
    A_1 = A + 25.0/216*k1 + 1408.0/2565*k3 + 2197.0/4104*k4 - 1.0/5*k5;

    C_  = C + 16.0/135*m1 + 6656.0/12825*m3 + 28561.0/56430*m4 - 9.0/50*m5 + 2.0/55*m6;
    C_1 = C + 25.0/216*m1 + 1408.0/2565*m3 + 2197.0/4104*m4 - 1.0/5*m5;

    vec R;
    double maxR;
    double rTmp;

    R = 1.0/dt*abs(A_ - A_1);
    maxR = max(R);
    for(uint i=0; i<C.n_cols; i++){
        R = 1.0/dt*abs(C_.col(i) - C_1.col(i));
        rTmp = max(R);
        if(rTmp > maxR)
            maxR = rTmp;
    }

    double delta = 0.84*pow(epsilon/maxR, 0.25);

    if(maxR <= epsilon){
        C = C_1;
        A = A_1;
        t += dt;
        dt = delta*dt;

        // Normalizing
        renormalize(C);
        A = A/sqrt(cdot(A,A));

        accepted = true;
    }else
    {
        dt = delta*dt;
        step--;
    }
    return accepted;
}
//------------------------------------------------------------------------------
