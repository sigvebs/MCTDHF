MCTDHF
================
An implementation of the Multiconfigurational Time Dependent Hartree-Fock method(MCTDHF) for solving time-dependent many-body problems in quantum mechanics. To compile the following libraries are needed:
- Armadillo
- Libconfig
- Lapack
- Blas
- FFTW

To compile run "qmake CONFIG+=default' and then 'make'. To run the program execute 'MCTDHF PATH-CONFIG'
where 'PATH-CONFIG' is the path to a configuration file.

An example config file:
```
#-------------------------------
# Example of configuration file
# for the MCTDHF program
#-------------------------------
systemSettings:
{
    version = "---";
    cleanFiles = true;
    plotResult = true;
    doTimeIntegration = false;

    saveToFileInterval = 100;
    filePath = "../DATA/";
};

system:
{
    shells      = 2;
    nParticles  = 2;
    dim         = 1;

    # 0 = Cartesian
    coordinateType = 0;

    conserveSpin    = false;
    spinValue       = 0;
};

spatialDiscretization:
{
    latticeRange = 10.0;
    nGrid = 8;

    # Differential operator:
    # 0 = "Finite Difference 1d"
    # 1 = "Finite Difference Five Point 1d"
    # 2 = "Spectral Method 1d"
    differentialOperator = 2;
};

ComplexTimeIntegration:
{
    # 0 = "Crank-Nicolson"
    # 1 = "Runge-Kutta 4"
    # 2 = "Runge-Kutta-Fehlberg"
    integrator = 1;
    dt = 0.005;
    N = 5000;

    rungeKuttaFehlberg:
    {
        epsilon = 0.00001;
    };
};

timeIntegration:
{
    # 0 = "Runge-Kutta 4"
    # 1 = "Runge-Kutta-Fehlberg"
    integrator = 1;
    dt = 0.01;
    N = 500;

    rungeKuttaFehlberg:
    {
        epsilon = 0.00001;
    };
};

interactionPotential:
{
    interactionType = 1;

    # "0"       "-epsilon|x-y|^2"
    interactionPotential:
    {
        epsilon = 0.2;
    };

    # "1"       " lambda/sqrt((x-y)^2 + a^2)"
    shieldedCoulombInteraction:
    {
        lambda = 1.0;
        a = 0.25; # Shielding parameter
    };
};

meanFieldIntegrator:
{
    # 0 = "Trapezodial"
    # 1 = "Low rank approximation"
    integratorType = 0;

    lowRankApproximation:
    {
        constEnd = 1.;
        constValue = 0.35;
        endValue = 1.0;
        epsilon = 1.0;
    };
};

wavefunction:
{
    # 0 = "Harmonic Oscillator"
    # 1 = "Hydrogen Like"
    basisType = 0;
};

oneBodyPotential:
{
    potential = [0];
    timeDepPotential = [2];

    # "0"       " 0.5*w^2x^2 "
    harmonicOscillatorBinding:
    {
         w = 0.25; # Strength of the confining potential.
    };

    # "1"       " Z/sqrt((x-y)^2 + b^2) "
    coulombInteractionNucleus:
    {
        Z = 2.0; # charge of the nucleus
        b = 0.7408; # Shielding paramter
    };

    # "2"       " x e0 sin(w t) "
    simpleLaser:
    {
        w = 8.0;  # Frequency of laser - measured in w from "harmonicOscillatorBinding"
        e0 = 1.0;  # Amplitude of laser
    };

    # "3"       " 1/(2d^2)(x - 0.5d)^2(x + 0.5d)^2"
    anharmonicDoubleWell:
    {
        d = 8.0;
    };
};
```
