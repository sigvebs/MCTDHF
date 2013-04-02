#ifndef MFTRAPEZOIDAL_H
#define MFTRAPEZOIDAL_H

// Local includes
#include <src/Interaction/MeanFieldIntegrator/meanfieldintegrator.h>

class MfTrapezoidal: public MeanFieldIntegrator
{
public:
    MfTrapezoidal(Config* cfg);
    virtual void initialize();
    virtual void integrate(const int q, const int r, const cx_mat &C, cx_vec &V2);
    virtual cx_double integrate(const int p, const int q, const int r, const int s, const cx_mat &C);
protected:
    mat Vxy;
    const double *Vxy_;
};

#endif // MFTRAPEZOIDAL_H
