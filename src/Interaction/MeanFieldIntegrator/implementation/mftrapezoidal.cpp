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
}
//------------------------------------------------------------------------------
cx_vec MfTrapezoidal::integrate(const int q, const int r, const cx_mat &C)
{
    Vtmp.zeros(nGrid);
    cx_double integral ;

    // Calulating the mean field V^qr
    for(int i=0; i<nGrid; i++){

        // Integrations using the trapezodial rule.
        integral = 0;
        for(int j=1; j<nGrid-1; j++){
            integral += conj(C(j,q))*Vxy(j,i)*C(j,r);
        }
        integral *=2;

        // Enpoints
        integral += conj(C(0,q))*Vxy(0,i)*C(0,r)
                + conj(C(nGrid-1,q))*Vxy(nGrid-1,i)*C(nGrid-1,r);

//        integral = cx_double(real(integral),0);
        Vtmp(i) = 0.5*integral;
    }
    return Vtmp;
}
//------------------------------------------------------------------------------
cx_double MfTrapezoidal::integrate(const int p, const int q, const int r, const int s, const cx_mat &C)
{
    // Integrations using the trapezodial rule.
    cx_double integral = 0;
    for(int i=1; i<nGrid-1; i++){
        integral += conj(C(i,p))*V2(q,s)(i)*C(i,r);
    }
    integral *= 2;

    // Enpoints
    integral += conj(C(0,p))*V2(q,s)(0)*C(0,r) + conj(C(nGrid-1,p))*V2(q,s)(nGrid-1)*C(nGrid-1,r);

    return 0.5*integral;
}
//------------------------------------------------------------------------------
