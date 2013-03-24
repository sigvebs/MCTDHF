#include "mftrapezoidal.h"

//------------------------------------------------------------------------------
MfTrapezoidal::MfTrapezoidal(Config *cfg):
    MeanFieldIntegrator(cfg)
{
}
//------------------------------------------------------------------------------
void MfTrapezoidal::initialize()
{
    Vxy = zeros(nGrid, nGrid);

    for(uint i=0; i<potential.size(); i++){
        Vxy += potential[i]->computeInteractionSpace();
    }
    Vtmp = zeros<cx_vec>(nGrid);
}
//------------------------------------------------------------------------------
cx_vec MfTrapezoidal::integrate(const int q, const int r, const cx_mat &C)
{
    cx_double integral;
    Vtmp.zeros();

    // Calulating the mean field V^qr
    for(int i=0; i<nGrid; i++){

        // Integrations using the trapezodial rule.
        // Enpoints
        integral = 0.5*(conj(C(0,q))*Vxy(0,i)*C(0,r)
                        + conj(C(nGrid-1,q))*Vxy(nGrid-1,i)*C(nGrid-1,r));

        for(int j=1; j<nGrid-1; j++){
            integral += conj(C(j,q))*Vxy(j,i)*C(j,r);
        }

        Vtmp(i) = integral;
    }
    return Vtmp;
}
//------------------------------------------------------------------------------
cx_double MfTrapezoidal::integrate(const int p, const int q, const int r, const int s, const cx_mat &C)
{
    // Integrations using the trapezodial rule.
    // Enpoints
    cx_double integral = 0.5*(conj(C(0,p))*V2(q,s)(0)*C(0,r) + conj(C(nGrid-1,p))*V2(q,s)(nGrid-1)*C(nGrid-1,r));

    for(int i=1; i<nGrid-1; i++){
        integral += conj(C(i,p))*V2(q,s)(i)*C(i,r);
    }

    return integral;
}
//------------------------------------------------------------------------------
