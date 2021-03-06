#include "mftrapezoidal.h"

//------------------------------------------------------------------------------
MfTrapezoidal::MfTrapezoidal(Config *cfg, const Grid &grid):
    MeanFieldIntegrator(cfg, grid)
{
}
//------------------------------------------------------------------------------
void MfTrapezoidal::initialize()
{
    Vxy = zeros(nGrid, nGrid);

    for(uint p=0; p<potential.size(); p++){
        for(int i=0; i<nGrid; i++){
            for(int j=0; j<nGrid; j++){
                Vxy(i, j) += potential[p]->evaluate(i, j);
            }
        }
    }

    cleanUp();
}
//------------------------------------------------------------------------------
void MfTrapezoidal::integrate(const int q, const int r, const cx_mat &C, cx_vec &V2)
{
    // Calulating the mean field V^qr
    cx_double integral;

#if 1 // Optimized version
    uint rnGrid = r*nGrid;
    uint qnGrid = q*nGrid;

    const cx_double *C_ = C.memptr();
    cx_double *V2_ = V2.memptr();

    // Calulating the mean field V^qr
    for(int i=0; i<nGrid; i++){
        uint inGrid = i*nGrid;

        // Integration using the trapezodial rule.
        integral = 0.5*(std::conj(C_[qnGrid])*Vxy[inGrid]*C_[rnGrid]
                        + std::conj( C_[nGrid-1 + q*nGrid])*Vxy[nGrid-1 + inGrid]*C_[nGrid-1 + rnGrid]);

        for(int j=1; j<nGrid-1; j++){
            integral += std::conj( C_[j + qnGrid] ) * Vxy[j + inGrid] * C_[j + rnGrid];
        }

        V2_[i] = integral;
    }
#else // Non-optimized version
    for(int i=0; i<nGrid; i++){

        // Integrations using the trapezodial rule.
        // Enpoints
        integral = 0.5*(std::conj(C(0,q))*Vxy(0,i)*C(0,r)
                        + std::conj(C(nGrid-1,q))*Vxy(nGrid-1,i)*C(nGrid-1,r));

        for(int j=1; j<nGrid-1; j++){
            integral += std::conj(C(j,q))*Vxy(j,i)*C(j,r);
        }

        V2(i) = integral;
    }
#endif
}
//------------------------------------------------------------------------------
cx_double MfTrapezoidal::integrate(const int p, const int q, const int r, const int s, const cx_mat &C)
{
    // Integrations using the trapezodial rule.

    // Enpoints
    cx_double integral = 0.5*(std::conj(C(0,p))*V2(q,s)(0)*C(0,r) + std::conj(C(nGrid-1,p))*V2(q,s)(nGrid-1)*C(nGrid-1,r));

    for(int i=1; i<nGrid-1; i++){
        integral += std::conj(C(i,p))*V2(q,s)(i)*C(i,r);
    }

    return integral;
}
//------------------------------------------------------------------------------
