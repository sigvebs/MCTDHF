#ifndef MCTDHFAPPLICATION_H
#define MCTDHFAPPLICATION_H

// Local includes
#include <src/includes/defines.h>

#include <src/Basis/basis.h>
#include <src/Basis/implementation/basisharmonicoscillator.h>
#include <src/Basis/implementation/basishydrogenlike.h>

#include <src/SlaterDeterminants/slaterdeterminants.h>

// Operators
#include <src/Interaction/interaction.h>
#include <src/OneParticleOperator/oneparticleoperator.h>

#include <src/OneParticleOperator/DifferentialOperator/differentialoperator.h>
#include <src/OneParticleOperator/DifferentialOperator/implementation/finitedifference1d.h>
#include <src/OneParticleOperator/DifferentialOperator/implementation/spectral1d.h>

#include <src/SlaterEquation/slaterequation.h>
#include <src/OrbitalEquation/orbitalequation.h>

#include <src/ComplexTimePropagation/complextimepropagation.h>

#include <src/ComplexTimePropagation/ComplexTimeIntegrator/complextimeintegrator.h>
#include <src/ComplexTimePropagation/ComplexTimeIntegrator/implementation/complextimecranknicholson.h>
#include <src/ComplexTimePropagation/ComplexTimeIntegrator/implementation/complextimerungekutta4.h>

// Libary incldues
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <libconfig.h++>
#include <sstream>

using namespace std;
using namespace libconfig;

class MctdhfApplication
{
public:
    MctdhfApplication(int *argc, char ***argv, string configFilename = "../config.cfg");
    void run();
protected:
    ComplexTimeIntegrator* setComplexTimeIntegrator();
    DifferentialOperator* setDifferentialOpertor();
    Basis* setBasis();
    Config cfg;
};

#endif // MCTDHFAPPLICATION_H
