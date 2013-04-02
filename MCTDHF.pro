CONFIG+=MPI
DEFINES += USE_MPI

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

cluster{
    INCLUDEPATH += /home/sigve/usr/local/include
    INCLUDEPATH += /home/sigve/usr/include
    LIBS += -L/home/sigve/usr/local/lib -lconfig++ -llapack -lblas -larmadillo -L/home/sigve/usr/lib -lfftw3 -lm
}

abel{
    QMAKE_CXX = mpicxx
    QMAKE_CXX_RELEASE = $$QMAKE_CXX
    QMAKE_CXX_DEBUG = $$QMAKE_CXX
    QMAKE_LINK = $$QMAKE_CXX
    QMAKE_CC = mpicc

    LIBS += -L/usit/abel/u1/sigve/usr/lib -lconfig++ -L/usit/abel/u1/sigve/usr/lib64 -larmadillo -L/usit/abel/u1/sigve/usr/lib -lfftw3 \
        -L$(MKLROOT)/lib/intel64 -lpthread  -mkl_intel_lp64 -liomp5

    QMAKE_CXXFLAGS += -DMKL_LP64 -I/cluster/software/VERSIONS/intel-2013.2/mkl/include -I$HOME/libs/armadillo/usr/include -openmp
    QMAKE_CXXFLAGS -= -m64
    INCLUDEPATH += /usit/abel/u1/sigve/usr/lib
    INCLUDEPATH += /usit/abel/u1/sigve/usr/lib64
    INCLUDEPATH += /usit/abel/u1/sigve/usr/include
}

MPI{
    QMAKE_CXX = mpicxx
    QMAKE_CXX_RELEASE = $$QMAKE_CXX
    QMAKE_CXX_DEBUG = $$QMAKE_CXX
    QMAKE_LINK = $$QMAKE_CXX
    QMAKE_CC = mpicc

    QMAKE_CFLAGS = $$system(mpicc --showme:compile)
    QMAKE_LFLAGS = $$system(mpicxx --showme:link)
    QMAKE_CXXFLAGS = $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
    QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
}

default{
    LIBS += -lconfig++ -llapack -lblas -larmadillo -lfftw3 -lm
}

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

QMAKE_CXXFLAGS += -std=c++0x

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
    src/OneParticleOperator/DifferentialOperator/implementation/finitedifference2d.cpp \
    src/OneParticleOperator/DifferentialOperator/implementation/finitedifferencefivepoint2d.cpp

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
    src/OneParticleOperator/DifferentialOperator/implementation/finitedifference2d.h \
    src/OneParticleOperator/DifferentialOperator/implementation/finitedifferencefivepoint2d.h

OTHER_FILES += \
    ../config.cfg \
    README.md \
    ../config.cfg
