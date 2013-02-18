#ifndef COMPLEXTIMERUNGEKUTTA4_H
#define COMPLEXTIMERUNGEKUTTA4_H

// Local includes
#include <src/ComplexTimePropagation/ComplexTimeIntegrator/complextimeintegrator.h>

class ComplexTimeRungeKutta4: public ComplexTimeIntegrator
{
public:
    ComplexTimeRungeKutta4(Config* cfg);
    virtual void stepForward();
protected:
    cx_mat renomralize(cx_mat C);
    cx_vec k1;
    cx_vec k2;
    cx_vec k3;
    cx_vec k4;

    cx_mat m1;
    cx_mat m2;
    cx_mat m3;
    cx_mat m4;

    double Eprev;
    double E;

    // TMP
    vec EVec;
    int step;
    cx_double i;
};

#endif // COMPLEXTIMERUNGEKUTTA4_H
