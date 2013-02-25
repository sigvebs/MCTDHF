#ifndef BASISHYDROGENLIKE_H
#define BASISHYDROGENLIKE_H

// Local includes
#include <src/Basis/basis.h>
#include <src/WaveFunction/wavefunction.h>
#include <src/WaveFunction/Implementations/hydrogenlike.h>

#include <src/includes/lib.h>

class BasisHydrogenLike: public Basis
{
public:
    BasisHydrogenLike(Config* cfg);
    virtual void createInitalDiscretization();
private:
    void discretization1d();

    double dx;
    mat x;
};

#endif // BASISHYDROGENLIKE_H
