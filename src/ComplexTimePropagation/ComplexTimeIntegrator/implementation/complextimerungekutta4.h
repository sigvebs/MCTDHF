#ifndef COMPLEXTIMERUNGEKUTTA4_H
#define COMPLEXTIMERUNGEKUTTA4_H

// Local includes
#include <src/ComplexTimePropagation/complextimepropagation.h>

class ComplexTimeRungeKutta4: public ComplexTimePropagation
{
public:
    ComplexTimeRungeKutta4(Config* cfg);
    virtual bool stepForward();
protected:
    cx_mat renormalize(cx_mat C);
    cx_vec k1;
    cx_vec k2;
    cx_vec k3;
    cx_vec k4;

    cx_mat m1;
    cx_mat m2;
    cx_mat m3;
    cx_mat m4;
};

#endif // COMPLEXTIMERUNGEKUTTA4_H
