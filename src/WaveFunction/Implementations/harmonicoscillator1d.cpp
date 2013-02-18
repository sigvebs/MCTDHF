#include "harmonicoscillator1d.h"

//------------------------------------------------------------------------------
HarmonicOscillator1d::HarmonicOscillator1d(Config *cfg, vec quantumNumbers):
    HarmonicOscillator(cfg, quantumNumbers)
{
    w = cfg->lookup("potential.w");
    n = quantumNumbers(1);
    coefficient = 1.0/(sqrt(pow(2,n)*factorial(n)))*pow(w/PI, 1.0/4);
}
//------------------------------------------------------------------------------
mat HarmonicOscillator1d::evaluate(const mat &x)
{
    double I;
    mat psi(x.n_rows, 1);
    for(int i=0;i<x.n_rows; i++){
        psi(i,0) = exp(-0.5*w*(x(i)) * x(i)) * hermitePolynomial(n, sqrtW * x(i));
    }
    return coefficient*psi;
}
//------------------------------------------------------------------------------
double HarmonicOscillator1d::getEnergy()
{
    return w*(n + 0.5);
}
//------------------------------------------------------------------------------
