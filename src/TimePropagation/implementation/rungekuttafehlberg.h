#ifndef RUNGEKUTTAFEHLBERG_H
#define RUNGEKUTTAFEHLBERG_H

// Local includes
#include <src/TimePropagation/timepropagation.h>

class RungeKuttaFehlberg: public TimePropagation
{
public:
    RungeKuttaFehlberg(Config* cfg);
    virtual bool stepForward();
protected:
    double epsilon;
};

#endif // RUNGEKUTTAFEHLBERG_H
