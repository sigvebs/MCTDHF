#ifndef ANHARMONICDOUBLEWELL_H
#define ANHARMONICDOUBLEWELL_H

// Local includes
#include <src/Potential/potential.h>

class AnharmonicDoubleWell: public Potential
{
public:
    AnharmonicDoubleWell(Config *cfg);
    virtual cx_vec evaluate(const cx_vec &psi, double t);
};

#endif // ANHARMONICDOUBLEWELL_H
