#include "mflowrankapproximation.h"

//------------------------------------------------------------------------------
MfLowRankApproximation::MfLowRankApproximation(Config *cfg, const Grid &grid):
    MeanFieldIntegrator(cfg, grid)
{
}
//------------------------------------------------------------------------------
void MfLowRankApproximation::initialize()
{
    double dx;
    double constValue;
    double endValue;
    double constEnd;
    double epsilon;


    dx = grid.DX;
    try {
        nOrbitals = cfg->lookup("spatialDiscretization.nSpatialOrbitals");
        constEnd = cfg->lookup("meanFieldIntegrator.lowRankApproximation.constEnd");
        constValue = cfg->lookup("meanFieldIntegrator.lowRankApproximation.constValue");
        endValue = cfg->lookup("meanFieldIntegrator.lowRankApproximation.endValue");
        epsilon = cfg->lookup("meanFieldIntegrator.lowRankApproximation.epsilon");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "MfLowRankApproximation::Error reading entry from config object." << endl;
        exit(EXIT_FAILURE);
    }

    int nConst = constEnd/dx;
    int constCenter = nGrid/2;
    Vxy = zeros(nGrid, nGrid);

    for(uint p=0; p<potential.size(); p++){
        for(int i=0; i<nGrid; i++){
            for(int j=0; j<nGrid; j++){
                Vxy(i, j) += potential[p]->evaluate(i, j);
            }
        }
    }

    // Using a simple discretization equal to the discretization of
    // the system.
    mat h = hExactSpatial();
//    mat h = hPiecewiseLinear();
    mat Q = eye(nGrid,nGrid);

    // Using a consant weight in the center of the potential
    // and a linear decrease from the center.
    vec g = gLinear(nGrid, constCenter, nConst, constValue, endValue);
    mat C = cMatrix(g, h, dx);

    vec lambda;
    mat eigvec;
    eig_sym(lambda, eigvec, inv(C.t())*Vxy*inv(C));
    mat Ut = C.t()*eigvec;

    mat QU = inv(Q)*Ut;

    // Sorting the eigenvalues by absoulte value and finding the number
    // of eigenvalues with abs(eigenval(i)) > epsilon
    uvec indices = sort_index(abs(lambda), 1);
    M = -1;
    for(uint m=0; m <lambda.n_rows; m++){
        cout << abs(lambda(indices(m))) << endl;
        if(abs(lambda(indices(m))) < epsilon){
            M = m;
            break;
        }
    }
//    cout << min(abs(lambda)) << endl;
//    cout << "hei" << endl;
//    cout << lambda << endl;
//    M = 63;
    if(M < 0){
        cerr << "MfLowRankApproximation:: no eigenvalues < epsilon found."
             << " Try setting epsilon to a higher number" << endl;
        exit(EXIT_FAILURE);
    }
    eigenval = zeros(M);
    for(int m=0; m <M; m++){
        eigenval(m) = lambda(indices(m));
    }

    // Calculating  the U matrix
    U = zeros(nGrid, M);
    for(int m=0; m < M; m++){
        for(uint j=0; j<h.n_rows; j++){
            U(j,m) = 0;
            for(uint i=0; i<h.n_cols; i++){
                U(j,m) += h(j,i)*QU(i,indices(m));
            }
        }
    }
    cout << "MfLowRankApproximation:: Trunction of eigenvalues at M = " << M << endl;

#if 1 // For testing the low rank approximation's accuracy
    mat appV = zeros(nGrid, nGrid);
    for(int i=0; i<nGrid; i++){
        for(int j=0; j<nGrid; j++){
            appV(i,j) = 0;
            for(int m=0; m<M; m++){
                appV(i,j) += eigenval(m)*U(i,m)*U(j,m);
            }
        }
    }
    mat diffV = abs(Vxy - appV);
    cout << "max_err = " << max(max(abs(Vxy - appV))) << endl;

    diffV.save("../DATA/diffV.mat");
    appV.save("../DATA/Vapp.mat");
    Vxy.save("../DATA/Vex.mat");
    cout << nGrid << endl;
//    exit(EXIT_SUCCESS);

#endif
    cout << "test" << endl;
    Vm = zeros<cx_vec>(M);
    Vqr = zeros<cx_vec>(nGrid);
}
//------------------------------------------------------------------------------
cx_vec MfLowRankApproximation::integrate(const int q, const int r, const cx_mat &C)
{
    Vqr.zeros();

    for(int m=0; m<M; m++){
        for(int j=0; j<nGrid; j++){
            Vm(m) += std::conj(C(j,q))*U(j,m)*C(j,r);
        }
    }

    for(int i=0; i<nGrid; i++){
        for(int m=0; m<M; m++){
            Vqr(i) += eigenval(m)*U(i,m)*Vm(m);
        }
    }

    return Vqr;
}
//------------------------------------------------------------------------------
void MfLowRankApproximation::integrate(const int q, const int r, const cx_mat &C, cx_vec &V2)
{
    Vm.zeros();
    V2.zeros();

    for(int m=0; m<M; m++){
        for(int j=0; j<nGrid; j++){
            Vm(m) += std::conj(C(j,q))*U(j,m)*C(j,r);
        }
    }

    for(int i=0; i<nGrid; i++){
        for(int m=0; m<M; m++){
            V2(i) += eigenval(m)*U(i,m)*Vm(m);
        }
    }
}
//------------------------------------------------------------------------------
cx_double MfLowRankApproximation::integrate(const int p, const int q, const int r, const int s, const cx_mat &C)
{
    // Integrations using the trapezodial rule.
    cx_double integral = 0;
    for(int i=1; i<nGrid-1; i++){
        integral += std::conj(C(i,p))*V2(q,s)(i)*C(i,r);
    }
    integral *= 2;

    // Enpoints
    integral += std::conj(C(0,p))*V2(q,s)(0)*C(0,r) + std::conj(C(nGrid-1,p))*V2(q,s)(nGrid-1)*C(nGrid-1,r);

    return 0.5*integral;
}
//------------------------------------------------------------------------------
mat MfLowRankApproximation::hExactSpatial()
{
    return eye(nGrid, nGrid);
}
//------------------------------------------------------------------------------
mat MfLowRankApproximation::hPiecewiseLinear()
{
//    mat h = zeros(ngrid, nGrid);
//    h.save("../DATA/h.mat", arma_ascii);
//    exit(EXIT_SUCCESS);
    return eye(nGrid, nGrid);
}
//------------------------------------------------------------------------------
vec MfLowRankApproximation::gLinear(int n, int constCenter,
                                    int nConst, double constValue,
                                    double endValue)
{
    vec weight = zeros(n);
    vec y = linspace(endValue, constValue , constCenter - nConst+1);

    int j = n-1;
    for(int i=0; i<constCenter - nConst; i++){
        weight(i) = weight(j) = y(i);
        j--;
    }
    for(int i=constCenter - nConst; i<=constCenter + nConst; i++){
        weight(i) = constValue;
    }
    return weight;
}
//------------------------------------------------------------------------------
mat MfLowRankApproximation::cMatrix(const vec &g, const mat &h, double dx)
{
    int n = h.n_rows;
    int nSpatial = g.n_rows;

    mat S = zeros(n, n);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            S(i,j) = 0;
            for(int m=0; m<nSpatial; m++){
                S(i,j) += g(m)*h(m,i)*h(m,j);
            }
        }
    }
    return chol(S);;
}
//------------------------------------------------------------------------------
