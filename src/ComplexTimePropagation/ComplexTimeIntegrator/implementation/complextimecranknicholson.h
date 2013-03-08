#ifndef COMPLEXTIMECRANKNICHOLSON_H
#define COMPLEXTIMECRANKNICHOLSON_H

// Local includes
//#include <src/ComplexTimePropagation/ComplexTimeIntegrator/complextimeintegrator.h>
#include <src/ComplexTimePropagation/complextimepropagation.h>

class ComplexTimeCrankNicholson: public ComplexTimePropagation
{
public:
    ComplexTimeCrankNicholson(Config* cfg);
    virtual bool stepForward();
protected:
    cx_mat H1;
    cx_mat H2;
    cx_mat I;
    int n;
};

#endif // COMPLEXTIMECRANKNICHOLSON_H
