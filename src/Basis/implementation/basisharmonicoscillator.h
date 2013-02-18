#ifndef BASISHARMONICOSCILLATOR_H
#define BASISHARMONICOSCILLATOR_H

// Local includes
#include <src/Basis/basis.h>
#include <src/WaveFunction/wavefunction.h>
#include <src/WaveFunction/Implementations/harmonicoscillator1d.h>

#include <src/includes/lib.h>

class BasisHarmonicOscillator: public Basis
{
public:
    BasisHarmonicOscillator(Config* cfg);
    virtual void createInitalDiscretization();
private:
    void discretization1d();

    double dx;
    mat x;
};

#endif // BASISHARMONICOSCILLATOR_H
