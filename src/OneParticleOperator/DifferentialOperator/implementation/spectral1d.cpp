#include "spectral1d.h"
//------------------------------------------------------------------------------
Spectral1d::Spectral1d(Config* cfg): DifferentialOperator(cfg)
{
}
//------------------------------------------------------------------------------
cx_vec Spectral1d::secondDerivative(const cx_vec &phi)
{
    cx_vec diff = phi;

//    int N = 2*nGrid;
//    double *input;
//    input=new double [N];
//    std::ifstream file("example.dat");
//    for(int i=0;i<N;i++) file >> input[i] ;

//    file.close();

//    fftw_complex *out;
//    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2)+1);

//    fftw_plan p;

//    p = fftw_plan_dft_r2c_1d(N,(double*)phi.memptr(),out,FFTW_ESTIMATE);

//    fftw_execute(p);

//    for(int i=0;i<(N/2)+1;i++)
//    {
////        cout << i <<"\t"<< out[i][0] << "\t" << out[i][1] <<std::endl;
//        diff(i) = cx_double(out[i][0], out[i][1]);
//    }

//    cout << diff << endl;
    cerr << "Spectral1d:: NOT YET IMPLMENTED" << endl;
    exit(0);
    return diff;
}
//------------------------------------------------------------------------------
