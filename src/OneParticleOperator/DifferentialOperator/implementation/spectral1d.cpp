#include "spectral1d.h"
//------------------------------------------------------------------------------
Spectral1d::Spectral1d(Config* cfg): DifferentialOperator(cfg)
{
}
//------------------------------------------------------------------------------
cx_vec Spectral1d::secondDerivative(const cx_vec &phi)
{
    cx_vec diff = phi;

//    fftw_complex *in, *out;
//     fftw_plan p;
//     ...
//     in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
//     out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
//     p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
//     .
//     fftw_execute(p); /* repeat as needed */

//     fftw_destroy_plan(p);
//     fftw_free(in); fftw_free(out);


    int N = nGrid;

    fftw_complex fx;
    memcpy( &fx, phi.memptr(), sizeof( fftw_complex ) );
    cout << "hei " << sizeof( fftw_complex ) << endl;
    exit(1);
    fftw_complex *out;
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    fftw_plan p;

//    p = fftw_plan_dft_1d(N, &fx, out,FFTW_FORWARD, FFTW_ESTIMATE);
//    p = fftw_plan_dft_r2c_1d(N,fx, out,FFTW_ESTIMATE);

//    fftw_execute(p);

    cout << phi << endl;
    for(int i=0;i<(N/2)+1;i++)
    {
//        cout << i <<"\t"<< out[i][0] << "\t" << out[i][1] <<std::endl;
        diff(i) = cx_double(out[i][0], out[i][1]);
    }

    cout << diff << endl;

    cerr << "Spectral1d:: NOT YET IMPLMENTED" << endl;
    diff.save("../DATA/fourier.mat", arma_ascii);
    exit(0);
    return diff;
}
//------------------------------------------------------------------------------
