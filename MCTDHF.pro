TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt


CONFIG += warn_on
CONFIG += wall

cluster{
    INCLUDEPATH += /home/sigve/usr/local/include
    INCLUDEPATH += /home/sigve/usr/include
}

SOURCES += main.cpp \
    src/mctdhfapplication.cpp \
    src/Basis/basis.cpp \
    src/SlaterDeterminants/slaterdeterminants.cpp \
    src/Interaction/interaction.cpp \
    src/WaveFunction/wavefunction.cpp \
    src/WaveFunction/Implementations/harmonicoscillator.cpp \
    src/includes/lib.cpp \
    src/OrbitalEquation/orbitalequation.cpp \
    src/SlaterEquation/slaterequation.cpp \
    src/includes/binaryoperations.cpp \
    src/Basis/implementation/basisharmonicoscillator.cpp \
    src/OneParticleOperator/oneparticleoperator.cpp \
    src/ComplexTimePropagation/complextimepropagation.cpp \
    src/ComplexTimePropagation/ComplexTimeIntegrator/implementation/complextimerungekutta4.cpp \
    src/ComplexTimePropagation/ComplexTimeIntegrator/implementation/complextimecranknicholson.cpp \
    src/OneParticleOperator/DifferentialOperator/differentialoperator.cpp \
    src/OneParticleOperator/DifferentialOperator/implementation/finitedifference1d.cpp \
    src/OneParticleOperator/DifferentialOperator/implementation/spectral1d.cpp \
    src/WaveFunction/Implementations/hydrogenlike.cpp \
    src/Basis/implementation/basishydrogenlike.cpp \
    src/Potential/potential.cpp \
    src/InteractionPotential/interactionpotential.cpp \
    src/InteractionPotential/implementation/harmonicoscillatorinteraction.cpp \
    src/InteractionPotential/implementation/screenedcoulombinteraction.cpp \
    src/OneParticleOperator/DifferentialOperator/implementation/finitedifferencefivepoint1d.cpp \
    src/Interaction/MeanFieldIntegrator/meanfieldintegrator.cpp \
    src/Interaction/MeanFieldIntegrator/implementation/mftrapezoidal.cpp \
    src/Interaction/MeanFieldIntegrator/implementation/mflowrankapproximation.cpp \
    src/ComplexTimePropagation/ComplexTimeIntegrator/implementation/complextimerungekuttafehlberg.cpp \
    src/Potential/implementation/harmonicoscillatoronebody.cpp \
    src/Potential/implementation/coulombinteractionnucleus.cpp \
    src/TimePropagation/timepropagation.cpp \
    src/TimePropagation/implementation/rungekutta4.cpp \
    src/Potential/implementation/simplelaser.cpp \
    src/Potential/implementation/anharmonicdoublewell.cpp \
    src/Basis/implementation/randomunitarymatrix.cpp \
    src/TimePropagation/implementation/rungekuttafehlberg.cpp \
    src/Grid/grid.cpp \
    src/OneParticleOperator/DifferentialOperator/implementation/finitedifference2d.cpp

HEADERS += \
    src/mctdhfapplication.h \
    src/includes/defines.h \
    src/Basis/basis.h \
    src/SlaterDeterminants/slaterdeterminants.h \
    src/Interaction/interaction.h \
    src/WaveFunction/wavefunction.h \
    src/WaveFunction/Implementations/harmonicoscillator.h \
    src/includes/lib.h \
    src/OrbitalEquation/orbitalequation.h \
    src/SlaterEquation/slaterequation.h \
    src/includes/binaryoperations.h \
    src/Basis/implementation/basisharmonicoscillator.h \
    src/OneParticleOperator/oneparticleoperator.h \
    src/ComplexTimePropagation/complextimepropagation.h \
    src/ComplexTimePropagation/ComplexTimeIntegrator/implementation/complextimerungekutta4.h \
    src/ComplexTimePropagation/ComplexTimeIntegrator/implementation/complextimecranknicholson.h \
    src/OneParticleOperator/DifferentialOperator/differentialoperator.h \
    src/OneParticleOperator/DifferentialOperator/implementation/finitedifference1d.h \
    src/OneParticleOperator/DifferentialOperator/implementation/spectral1d.h \
    src/WaveFunction/Implementations/hydrogenlike.h \
    src/Basis/implementation/basishydrogenlike.h \
    src/Potential/potential.h \
    src/InteractionPotential/interactionpotential.h \
    src/InteractionPotential/implementation/harmonicoscillatorinteraction.h \
    src/InteractionPotential/implementation/screenedcoulombinteraction.h \
    src/OneParticleOperator/DifferentialOperator/implementation/finitedifferencefivepoint1d.h \
    src/Interaction/MeanFieldIntegrator/meanfieldintegrator.h \
    src/Interaction/MeanFieldIntegrator/implementation/mftrapezoidal.h \
    src/Interaction/MeanFieldIntegrator/implementation/mflowrankapproximation.h \
    src/ComplexTimePropagation/ComplexTimeIntegrator/implementation/complextimerungekuttafehlberg.h \
    src/Potential/implementation/harmonicoscillatoronebody.h \
    src/Potential/implementation/coulombinteractionnucleus.h \
    src/TimePropagation/timepropagation.h \
    src/TimePropagation/implementation/rungekutta4.h \
    src/Potential/implementation/simplelaser.h \
    src/Potential/implementation/anharmonicdoublewell.h \
    src/Basis/implementation/randomunitarymatrix.h \
    src/TimePropagation/implementation/rungekuttafehlberg.h \
    src/Grid/grid.h \
    src/OneParticleOperator/DifferentialOperator/implementation/finitedifference2d.h

OTHER_FILES += \
    ../config.cfg \
    README.md

default{
    LIBS += -lconfig++ -llapack -lblas -larmadillo -lfftw3 -lm
}

cluster{
    LIBS += -L/home/sigve/usr/local/lib -lconfig++ -llapack -lblas -larmadillo -L/home/sigve/usr/lib -lfftw3 -lm
}

home{
    LIBS += -lconfig++ -llapack -lblas -larmadillo -lfftw3 -lm
}
UIO_noIntel{
    LIBS += -fopenmp -lconfig++ -llapack -lblas -larmadillo -lfftw3 -lm
}
UIO{
    LIBS += -lconfig++ -llapack -lblas -larmadillo -lfftw3 -lm \
#    -L$(MKLROOT)/lib/intel64 -mkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread \
}

QMAKE_CXXFLAGS_DEBUG += -std=c++0x
QMAKE_CXXFLAGS_RELEASE += -std=c++0x

CONFIG(debug, debug|release) {
#    DEFINES += DEBUG
}

release {
    # Remoing other O flags
    QMAKE_CXXFLAGS_RELEASE -= -O
    QMAKE_CXXFLAGS_RELEASE -= -O1
    QMAKE_CXXFLAGS_RELEASE -= -O2

    # add the desired -O3 if not present
    QMAKE_CXXFLAGS_RELEASE *= -O3
    DEFINES += ARMA_NO_DEBUG
}
