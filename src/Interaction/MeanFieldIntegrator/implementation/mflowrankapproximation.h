#ifndef MFLOWRANKAPPROXIMATION_H
#define MFLOWRANKAPPROXIMATION_H

// Local includes
#include <src/Interaction/MeanFieldIntegrator/meanfieldintegrator.h>

class MfLowRankApproximation: public MeanFieldIntegrator
{
public:
    MfLowRankApproximation(Config* cfg);
    virtual void initialize();
    virtual cx_vec integrate(const int q, const int r, const cx_mat &C);
    virtual cx_double integrate(const int p, const int q, const int r, const int s, const cx_mat &C);
protected:
    mat U;
    mat V_;
    int M;
    vec eigenval;
    cx_vec Vtmp;

    // Local functions
    mat hExactSpatial();
    mat hPiecewiseLinear();
    vec gLinear(int n, int constCenter, int nConst, double constValue, double endValue);
    mat cMatrix(const vec &g, const mat &h, double dx);
};

#endif // MFLOWRANKAPPROXIMATION_H
