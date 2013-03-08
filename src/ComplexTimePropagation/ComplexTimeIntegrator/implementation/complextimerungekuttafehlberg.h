#ifndef COMPLEXTIMERUNGEKUTTAFEHLBERG_H
#define COMPLEXTIMERUNGEKUTTAFEHLBERG_H

// Local includes
#include <src/ComplexTimePropagation/complextimepropagation.h>

class ComplexTimeRungeKuttaFehlberg: public ComplexTimePropagation
{
public:
    ComplexTimeRungeKuttaFehlberg(Config* cfg);
    virtual bool stepForward();
protected:
    double epsilon;

    cx_mat renormalize(cx_mat C);
    cx_vec k1;
    cx_vec k2;
    cx_vec k3;
    cx_vec k4;
    cx_vec k5;
    cx_vec k6;

    cx_mat m1;
    cx_mat m2;
    cx_mat m3;
    cx_mat m4;
    cx_mat m5;
    cx_mat m6;

    cx_vec A_;
    cx_mat C_;
    cx_vec A_1;
    cx_mat C_1;
};

#endif // COMPLEXTIMERUNGEKUTTAFEHLBERG_H
