#ifndef RUNGEKUTTA4_H
#define RUNGEKUTTA4_H

// Local includes
#include <src/TimePropagation/timepropagation.h>

class RungeKutta4: public TimePropagation
{
public:
    RungeKutta4(Config* cfg);
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

#endif // RUNGEKUTTA4_H
