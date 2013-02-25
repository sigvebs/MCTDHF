#include "hydrogenlike.h"

//------------------------------------------------------------------------------
HydrogenLike::HydrogenLike(Config *cfg, vec quantumNumbers):
    Wavefunction(cfg, quantumNumbers)
{
    n = quantumNumbers(1);
}
//------------------------------------------------------------------------------
mat HydrogenLike::evaluate(const mat &x)
{
    mat psi(x.n_rows, 1);
    double r;
     for(uint i=0;i<x.n_rows; i++){
         r = sqrt(x(i)*x(i));
         psi(i,0) = laguerrePolynomial(n,x(i))*exp(-abs(x(i)));
//         psi(i,0) = laguerrePolynomial(m,r)*r*exp(-0.5*r);
     }
     return psi;
}
//------------------------------------------------------------------------------
double HydrogenLike::laguerrePolynomial(const int n, const double x)
{
    double pol;

    switch(n){
    case 0:
        pol = 1;
        break;
    case 1:
        pol = -x + 1;
        break;
    case 2:
        pol = 0.5*(pow(x,2) -4*x +2);
        break;
    case 3:
        pol = 1.0/6.0*(-pow(x,3) + 9*pow(x,2) - 18*x + 6);
        break;
    case 4:
        pol = 1.0/24.0*(pow(x,4) - 16*pow(x,3) + 72*pow(x,2) - 96*x + 24);
        break;
    case 5:
        pol = 1.0/120.0*(-pow(x,5) + 25*pow(x,4) - 200*pow(x,3) + 600*pow(x,2) + 120);
        break;
    case 6:
        pol = 1.0/720.0*(pow(x,6) - 36*pow(x,5) + 450*pow(x,4) - 2400*pow(x,3) + 5400*pow(x,2) - 4320*x +720);
        break;
    default:
        cerr << "HydrogenLike::laguerrePolynomial::Laguerre Polynomial of degree "
             << n << " has not been implmented" << endl;
        exit(EXIT_FAILURE);
        break;
    }
    return pol;
}
//------------------------------------------------------------------------------
