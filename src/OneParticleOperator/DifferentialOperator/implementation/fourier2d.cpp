#include "fourier2d.h"

//------------------------------------------------------------------------------
Fourier2d::Fourier2d(Config* cfg, const Grid &grid):
    DifferentialOperator(cfg, grid)
{
    try{
        nGridX = cfg->lookup("spatialDiscretization.twoD.nGridX");
        nGridY = cfg->lookup("spatialDiscretization.twoD.nGridY");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "DifferentialOperator::DifferentialOperator(Config *cfg)"
             << "::Error reading from config object." << endl;
        exit(EXIT_FAILURE);
    }

    // Setting up the frequencies.

    //------------------------------
    // X
    vec k_x = vec(nGridX);

    for(int i=0; i<nGridX/2; i++) {
        k_x[i] = i;
    }
    for(int i=nGridX/2; i<nGridX;i++) {
        k_x[i] = - (nGridX - i);
    }

    k_x *= 2*PI/(dx*(nGridX));
    k_x = -k_x%k_x; // Included the scaling 1.0/N here

    //------------------------------
    // Y
    vec k_y = vec(nGridY);

    for(int i=0; i<nGridY/2; i++) {
        k_y[i] = i;
    }
    for(int i=nGridY/2; i<nGridY;i++) {
        k_y[i] = - (nGridY - i);
    }

    k_y *= 2*PI/(dy*(nGridY));
    k_y = -k_y%k_y; // Included the scaling 1.0/N here
    //------------------------------

    // Creating the total wavespace
    k = vec(nGrid);
    int c = 0;
    for(int i=0; i<nGridX; i++) {
        for(int j=0; j<nGridY; j++) {
            k(c++) = k_x(i) + k_y(j);
        }
    }
    k /= nGrid;

    diff = cx_vec(nGrid);

    forward = fftw_plan_dft_2d(nGridX, nGridY, (fftw_complex*)diff.memptr(), (fftw_complex*)diff.memptr(), FFTW_FORWARD, FFTW_ESTIMATE);
    backward = fftw_plan_dft_2d(nGridX, nGridY, (fftw_complex*)diff.memptr(), (fftw_complex*)diff.memptr(), FFTW_BACKWARD, FFTW_ESTIMATE);
}
//------------------------------------------------------------------------------
cx_vec Fourier2d::secondDerivative(const cx_vec &phi)
{
    diff = phi;
    fftw_execute(forward);
    diff = k % diff;
    fftw_execute(backward);
    return diff;
}
//------------------------------------------------------------------------------
Fourier2d::~Fourier2d()
{
    fftw_destroy_plan(forward);
    fftw_destroy_plan(backward);
}
//------------------------------------------------------------------------------
