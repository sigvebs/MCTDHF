TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

#INCLUDEPATH += /home/sigve/usr/local/include

SOURCES += main.cpp \
    src/mctdhfapplication.cpp \
    src/Basis/basis.cpp \
    src/SlaterDeterminants/slaterdeterminants.cpp \
    src/Interaction/interaction.cpp \
    src/WaveFunction/wavefunction.cpp \
    src/WaveFunction/Implementations/harmonicoscillator1d.cpp \
    src/WaveFunction/Implementations/harmonicoscillator.cpp \
    src/includes/lib.cpp \
    src/OrbitalEquation/orbitalequation.cpp \
    src/SlaterEquation/slaterequation.cpp \
    src/includes/binaryoperations.cpp \
    src/Basis/implementation/basisharmonicoscillator.cpp \
    src/OneParticleOperator/oneparticleoperator.cpp \
    src/ComplexTimePropagation/complextimepropagation.cpp \
    src/ComplexTimePropagation/ComplexTimeIntegrator/complextimeintegrator.cpp \
    src/ComplexTimePropagation/ComplexTimeIntegrator/implementation/complextimerungekutta4.cpp \
    src/ComplexTimePropagation/ComplexTimeIntegrator/implementation/complextimecranknicholson.cpp \
    src/OneParticleOperator/DifferentialOperator/differentialoperator.cpp \
    src/OneParticleOperator/DifferentialOperator/implementation/finitedifference1d.cpp \
    src/OneParticleOperator/DifferentialOperator/implementation/spectral1d.cpp \
    src/WaveFunction/Implementations/hydrogenlike.cpp \
    src/Basis/implementation/basishydrogenlike.cpp \
    src/Interaction/implementation/screenedcoulomb.cpp \
    src/Potential/potential.cpp \
    src/InteractionPotential/interactionpotential.cpp \
    src/InteractionPotential/implementation/harmonicoscillatorinteraction.cpp \
    src/InteractionPotential/implementation/screenedcoulombinteraction.cpp \
    src/OneParticleOperator/DifferentialOperator/implementation/finitedifferencefivepoint1d.cpp \
    src/Interaction/MeanFieldIntegrator/meanfieldintegrator.cpp \
    src/Interaction/MeanFieldIntegrator/implementation/mftrapezoidal.cpp \
    src/Interaction/MeanFieldIntegrator/implementation/mflowrankapproximation.cpp

HEADERS += \
    src/mctdhfapplication.h \
    src/includes/defines.h \
    src/Basis/basis.h \
    src/SlaterDeterminants/slaterdeterminants.h \
    src/Interaction/interaction.h \
    src/WaveFunction/wavefunction.h \
    src/WaveFunction/Implementations/harmonicoscillator1d.h \
    src/WaveFunction/Implementations/harmonicoscillator.h \
    src/includes/lib.h \
    src/OrbitalEquation/orbitalequation.h \
    src/SlaterEquation/slaterequation.h \
    src/includes/binaryoperations.h \
    src/Basis/implementation/basisharmonicoscillator.h \
    src/OneParticleOperator/oneparticleoperator.h \
    src/ComplexTimePropagation/complextimepropagation.h \
    src/ComplexTimePropagation/ComplexTimeIntegrator/complextimeintegrator.h \
    src/ComplexTimePropagation/ComplexTimeIntegrator/implementation/complextimerungekutta4.h \
    src/ComplexTimePropagation/ComplexTimeIntegrator/implementation/complextimecranknicholson.h \
    src/OneParticleOperator/DifferentialOperator/differentialoperator.h \
    src/OneParticleOperator/DifferentialOperator/implementation/finitedifference1d.h \
    src/OneParticleOperator/DifferentialOperator/implementation/spectral1d.h \
    src/WaveFunction/Implementations/hydrogenlike.h \
    src/Basis/implementation/basishydrogenlike.h \
    src/Interaction/implementation/screenedcoulomb.h \
    src/Potential/potential.h \
    src/InteractionPotential/interactionpotential.h \
    src/InteractionPotential/implementation/harmonicoscillatorinteraction.h \
    src/InteractionPotential/implementation/screenedcoulombinteraction.h \
    src/OneParticleOperator/DifferentialOperator/implementation/finitedifferencefivepoint1d.h \
    src/Interaction/MeanFieldIntegrator/meanfieldintegrator.h \
    src/Interaction/MeanFieldIntegrator/implementation/mftrapezoidal.h \
    src/Interaction/MeanFieldIntegrator/implementation/mflowrankapproximation.h

OTHER_FILES += \
    ../config.cfg

LIBS += -lconfig++ -larmadillo -llapack -lblas -lfftw3 -lm

QMAKE_CXXFLAGS_DEBUG += -std=c++0x
QMAKE_CXXFLAGS_RELEASE += -std=c++0x

CONFIG(debug, debug|release) {
    DEFINES += DEBUG
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

# For creating a version.txt file in the build directory
version.target = version
version.commands = python $$PWD/version.py $$PWD $$OUT_PWD
QMAKE_EXTRA_TARGETS += version
PRE_TARGETDEPS += version
