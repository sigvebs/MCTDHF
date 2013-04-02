#include "hydrogenlike.h"

//------------------------------------------------------------------------------
HydrogenLike::HydrogenLike(Config *cfg, vec quantumNumbers):
    Wavefunction(cfg, quantumNumbers)
{
}

//------------------------------------------------------------------------------
double HydrogenLike::evaluate(const vec &r)
{
    int dim = r.n_rows;
    double laguerre = 1;
    double rNorm = norm(r,2);

    for(int d=0; d<dim; d++){
        laguerre *= laguerrePolynomial(quantumNumbers(d+1), r(d));
    }

    return laguerre*exp(-rNorm);
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
