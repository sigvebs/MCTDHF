#ifndef COMPLEXTIMECRANKNICHOLSON_H
#define COMPLEXTIMECRANKNICHOLSON_H

// Local includes
#include <src/ComplexTimePropagation/ComplexTimeIntegrator/complextimeintegrator.h>

class ComplexTimeCrankNicholson: public ComplexTimeIntegrator
{
public:
    ComplexTimeCrankNicholson(Config* cfg);
    virtual void stepForward();
protected:
    cx_mat H1;
    cx_mat H2;
    cx_mat I;
    int n;
};

#endif // COMPLEXTIMECRANKNICHOLSON_H
