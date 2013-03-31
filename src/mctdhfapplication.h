#ifndef MCTDHFAPPLICATION_H
#define MCTDHFAPPLICATION_H

// Local includes
#include <src/includes/defines.h>

#include <src/Grid/grid.h>

#include <src/Basis/basis.h>
#include <src/Basis/implementation/basisharmonicoscillator.h>
#include <src/Basis/implementation/basishydrogenlike.h>
#include <src/Basis/implementation/randomunitarymatrix.h>

#include <src/SlaterDeterminants/slaterdeterminants.h>

// Operators
#include <src/Interaction/interaction.h>
#include <src/OneParticleOperator/oneparticleoperator.h>

// Interactions
#include <src/InteractionPotential/interactionpotential.h>
#include <src/InteractionPotential/implementation/harmonicoscillatorinteraction.h>
#include <src/InteractionPotential/implementation/screenedcoulombinteraction.h>

// Mean Field integration
#include <src/Interaction/MeanFieldIntegrator/meanfieldintegrator.h>
#include <src/Interaction/MeanFieldIntegrator/implementation/mftrapezoidal.h>
#include <src/Interaction/MeanFieldIntegrator/implementation/mflowrankapproximation.h>

// Differential Operators
#include <src/OneParticleOperator/DifferentialOperator/differentialoperator.h>
#include <src/OneParticleOperator/DifferentialOperator/implementation/finitedifference1d.h>
#include <src/OneParticleOperator/DifferentialOperator/implementation/finitedifference2d.h>
#include <src/OneParticleOperator/DifferentialOperator/implementation/finitedifferencefivepoint1d.h>
#include <src/OneParticleOperator/DifferentialOperator/implementation/spectral1d.h>

// One body operators
#include <src/Potential/potential.h>
#include <src/Potential/implementation/harmonicoscillatoronebody.h>
#include <src/Potential/implementation/coulombinteractionnucleus.h>
#include <src/Potential/implementation/simplelaser.h>
#include <src/Potential/implementation/anharmonicdoublewell.h>

#include <src/SlaterEquation/slaterequation.h>
#include <src/OrbitalEquation/orbitalequation.h>

// Imaginary time integration
#include <src/ComplexTimePropagation/complextimepropagation.h>

#include <src/ComplexTimePropagation/ComplexTimeIntegrator/implementation/complextimerungekutta4.h>
#include <src/ComplexTimePropagation/ComplexTimeIntegrator/implementation/complextimerungekuttafehlberg.h>

// Time integration
#include <src/TimePropagation/timepropagation.h>

#include <src/TimePropagation/implementation/rungekutta4.h>
#include <src/TimePropagation/implementation/rungekuttafehlberg.h>


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
    MctdhfApplication(string configFilename = "../config.cfg");
    void run();
protected:
    void setInteractionPotentials(Interaction &V, const Grid &grid);
    void setOneBodyPotentials(SingleParticleOperator &h, const Grid &grid);
    void setTimeDepOneBodyPotentials(SingleParticleOperator &h, const Grid &grid);
    TimePropagation* setTimeIntegrator();
    ComplexTimePropagation *setComplexTimeIntegrator();
    DifferentialOperator* setDifferentialOpertor(const Grid &grid);
    MeanFieldIntegrator* setMeanFieldIntegrator();
    Basis* setBasis();
    Config cfg;

    bool loadDataset;
    string loadDatasetPath;
    bool doTimeIntegration;
    bool doComplexTimeIntegration;
};

#endif // MCTDHFAPPLICATION_H
