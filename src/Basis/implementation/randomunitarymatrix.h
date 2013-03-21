#ifndef RANDOMUNITARYMATRIX_H
#define RANDOMUNITARYMATRIX_H

// Library includes
#include <stdlib.h>
#include <time.h>

// Local includes
#include <src/Basis/basis.h>
#include <src/WaveFunction/wavefunction.h>
#include <src/WaveFunction/Implementations/harmonicoscillator1d.h>

#include <src/includes/lib.h>

class RandomUnitaryMatrix: public Basis
{
public:
    RandomUnitaryMatrix(Config* cfg);
    virtual void createInitalDiscretization();
private:
    void discretization1d();
};

#endif // RANDOMUNITARYMATRIX_H
