#include "spectral1d.h"
//------------------------------------------------------------------------------
Spectral1d::Spectral1d(Config* cfg, const Grid &grid):
    DifferentialOperator(cfg, grid)
{
    // Setting up the frequencies.
    k = vec(nGrid);
    for(int i=0; i<nGrid/2; i++) {
        k[i] = i;
    }
    for(int i=nGrid/2; i<nGrid;i++) {
        k[i] = - (nGrid - i);
    }

    k *= 2*PI/(dx*(nGrid));
    k = -k%k/nGrid; // Included the scaling 1.0/N here

    diff = cx_vec(nGrid);

    forward = fftw_plan_dft_1d(nGrid, (fftw_complex*)diff.memptr(), (fftw_complex*)diff.memptr(), FFTW_FORWARD, FFTW_ESTIMATE);
    backward = fftw_plan_dft_1d(nGrid, (fftw_complex*)diff.memptr(), (fftw_complex*)diff.memptr(), FFTW_BACKWARD, FFTW_ESTIMATE);
}
//------------------------------------------------------------------------------
cx_vec Spectral1d::secondDerivative(const cx_vec &phi)
{
    diff = phi;

    fftw_execute(forward);
    diff = k % diff;
    fftw_execute(backward);

    return diff;
}
//------------------------------------------------------------------------------
Spectral1d::~Spectral1d()
{
    fftw_destroy_plan(forward);
    fftw_destroy_plan(backward);
}
//------------------------------------------------------------------------------
