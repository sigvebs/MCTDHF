#ifndef COMPLEXTIMERUNGEKUTTA4_H
#define COMPLEXTIMERUNGEKUTTA4_H

// Local includes
#include <src/ComplexTimePropagation/complextimepropagation.h>

class ComplexTimeRungeKutta4: public ComplexTimePropagation
{
public:
    ComplexTimeRungeKutta4(Config* cfg);
    virtual bool stepForward();
};

#endif // COMPLEXTIMERUNGEKUTTA4_H
