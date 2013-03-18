#ifndef MFTRAPEZOIDAL_H
#define MFTRAPEZOIDAL_H

// Local includes
#include <src/Interaction/MeanFieldIntegrator/meanfieldintegrator.h>

class MfTrapezoidal: public MeanFieldIntegrator
{
public:
    MfTrapezoidal(Config* cfg);
    virtual void initialize();
    virtual cx_vec integrate(const int q, const int r, const cx_mat &C);
    virtual cx_double integrate(const int p, const int q, const int r, const int s, const cx_mat &C);
protected:
    mat Vxy;
//    cx_vec Vtmp;
};

#endif // MFTRAPEZOIDAL_H
