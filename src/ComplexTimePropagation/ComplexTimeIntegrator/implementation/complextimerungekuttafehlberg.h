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

};

#endif // COMPLEXTIMERUNGEKUTTAFEHLBERG_H
