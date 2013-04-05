#include "orbitalequation.h"

//------------------------------------------------------------------------------
OrbitalEquation::OrbitalEquation(Config *cfg,
                                 vector<bitset<BITS> > slaterDeterminants,
                                 Interaction *V,
                                 SingleParticleOperator *h):
    cfg(cfg)
{
    this->slaterDeterminants = slaterDeterminants;
    this->V = V;
    this->h = h;

    try {
        nParticles = cfg->lookup("system.nParticles");
        nGrid = cfg->lookup("spatialDiscretization.nGrid");
        nOrbitals = cfg->lookup("spatialDiscretization.nSpatialOrbitals");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "OrbitalEquation::Error reading from config object." << endl;
        exit(EXIT_FAILURE);
    }
    nSlaterDeterminants = slaterDeterminants.size();
    invRho = cx_mat(nOrbitals, nOrbitals);

    U = zeros<cx_mat>(nGrid, nOrbitals);
    rightHandSide = zeros<cx_mat>(nGrid, nOrbitals);

//    Q = cx_mat(nGrid, nGrid);
    rho1 = zeros<cx_mat>(2*nOrbitals, 2*nOrbitals); // TMP

    myRank = 0;
    nNodes = 1;
#ifdef USE_MPI
    // MPI-------------------------------------------
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &nNodes);
#endif
    sizeRij = ivec(nNodes);
    int tot = nGrid*nOrbitals;

    allRH = imat(nGrid, nOrbitals);
    int node = 0;
    int s = 0;
    for (int j = 0; j < nOrbitals; j++) {
        for (int i = 0; i < nGrid; i++) {
            if (myRank == node){
                myRij.push_back(pair<int, int>(i,j));
            }
            allRH(i,j) = node;
            s++;
            if(s >= BLOCK_SIZE(node, nNodes, tot)){
                sizeRij(node) = BLOCK_SIZE(node, nNodes, tot);
                s = 0;
                node++;
            }
        }
    }
}
//------------------------------------------------------------------------------
const cx_mat &OrbitalEquation::computeRightHandSide(const cx_mat &C, const cx_vec &A)
{
    this->A = &A;
    hC = &(h->getHspatial());
    rho2.clear();

    computeOneParticleReducedDensity();
    invRho = inv(invRho);
    computeTwoParticleReducedDensity();
    computeUMatrix(C);

    //--------------------------------------------------------------------------
    // Hardcoding of the matrirx-matrix product: R.H = Q*U
    //--------------------------------------------------------------------------
    cx_double RH_ij;
    cx_double Qik;
    const cx_double *C_ = C.memptr();
    cx_double *U_ = U.memptr();
    cx_double *RH_ = rightHandSide.memptr();

    int i, j;
    for(pair<int,int> ij :myRij){
        i = ij.first;
        j = ij.second;

        RH_ij = 0;
        for(int k=0; k<nGrid; k++){

            // Calculating the projector
            Qik = 0;
            if(i == k)
                Qik = 1;
            for(int l=0; l<nOrbitals; l++){
                Qik -= C_[i + l*nGrid]*conj(C_[k + l*nGrid]);
            }
            // END projector
            RH_ij += Qik*U_[k + j*nGrid];
        }
        RH_[i + j*nGrid] = RH_ij;
    }

#ifdef USE_MPI
    i = 0;
    for(uint n=0; n<sizeRij.n_elem; n++){
        MPI_Bcast( &RH_[i], sizeRij(n) , MPI_DOUBLE_COMPLEX, n, MPI_COMM_WORLD );
        i += sizeRij(n);
    }
#endif
    return rightHandSide;
}
//------------------------------------------------------------------------------
void OrbitalEquation::computeUMatrix(const cx_mat &C)
{
    cx_vec Ui(nGrid);
    cx_vec rho_iqrs_Cs(nGrid);

    for(int j=0; j<nOrbitals; j++){
        U.col(j) = hC->col(j);

        for(int i=0; i<nOrbitals; i++){
            Ui.zeros();

            for(int q=0; q<nOrbitals; q++){
                for(int r=0; r<nOrbitals; r++){
                    rho_iqrs_Cs.zeros();

                    for(int s=0; s<nOrbitals; s++){
                        rho_iqrs_Cs += findRho2(i,q,r,s) * C.col(s);
                    }
                    Ui += V->meanField(q,r)%rho_iqrs_Cs;
                }
            }
            U.col(j) += invRho(j,i)*Ui;
        }
    }

#ifdef DEBUG
    cout << "OrbitalEquation::computeUMatrix()" << endl;
    // cout << "U = " << endl << U << endl;
    // exit(0);
    // Hardcoded the case for 2 electrons, 1 slater determinant, one spatial orbital.
    // rightHandSide = (I - C.col(0)*C.col(0).t())*(hC->col(0) + diagmat(V->meanField(0,0))*C.col(0));
#endif
}
//------------------------------------------------------------------------------
void OrbitalEquation::computeProjector(const cx_mat &C)
{
#if 1
    cx_double Qtmp;
    for(int m=0; m<nGrid; m++){
        for(int n=0; n<nGrid; n++){
            Qtmp = 0;
            for(int i=0; i<nOrbitals; i++){
                Qtmp -= C(m,i)*conj(C(n,i));
            }
            Q(m,n) = Qtmp;
        }
    }
#else
    // NEEDS TO BE UPDATED - not working
    // Slightly changing the projector to ensure orthonormality is conserved
    // Problems arise due to numerical inaccuracies.
    cx_mat O = zeros<cx_mat>(nSpatialOrbitals,nSpatialOrbitals);

    for(int i=0; i<nSpatialOrbitals; i++){
        for(int l=0; l<nSpatialOrbitals; l++){
            O(i,l) = cdot(C.col(i), C.col(l));
        }
    }
    O = inv(O);

    for(int i=0; i<nSpatialOrbitals; i++){
        for(int l=0; l<nSpatialOrbitals; l++){
            P += C.col(i)*C.col(l).t()*O(i,l);
        }
    }
    P = -P;
#endif

    // Adding the identity matrix
    for(int i=0; i<nGrid; i++){
        Q(i,i) += 1;
    }
#ifdef DEBUG
    cout << "OrbitalEquation::computeProjector()" << endl;
//    cout << "C0 * P* C0 = " << abs(C->col(0).t()*P*C->col(0)) << endl;
//    cout << "C0 * P * C1 = " << abs(C->col(0).t()*P*C->col(1)) << endl;
#endif
}
//------------------------------------------------------------------------------
void OrbitalEquation::computeOneParticleReducedDensity()
{
    // All possible spatial orbitals
    for(int i=0; i< nOrbitals; i++){
        for(int j=i; j<nOrbitals; j++){
            invRho(i,j)  = reducedOneParticleOperator(2*i,2*j);
            invRho(i,j)  += reducedOneParticleOperator(2*i+1,2*j+1);
            invRho(j,i) = conj(invRho(i,j));
        }
    }

#ifdef DEBUG
//#if 1
    cout << "OrbitalEquation::computeOneParticleReducedDensity()" << endl;
    cout << "Trace(rho) = " << real(trace(invRho)) << " nParticles = " << nParticles << endl;
#endif


#if 0
    //---------------------------
    // Regularization
    //---------------------------
    const double e = 1e-10;
    cx_mat rho_U;
    vec rho_lambda;
    eig_sym(rho_lambda, rho_U, invRho);

    rho_lambda = rho_lambda + e*exp(-rho_lambda*(1.0/e));
    cx_mat X = rho_U*diagmat(1.0/rho_lambda)*rho_U.t();

    invRho += e*X;
    //---------------------------
#endif

}
//------------------------------------------------------------------------------
void OrbitalEquation::computeTwoParticleReducedDensity()
{
    cx_double value;

    // All different
    for(int p=0; p<nOrbitals; p++){
        for(int q=0; q<nOrbitals; q++){
            for(int r=0; r<nOrbitals; r++){
                for(int s=0; s<nOrbitals; s++){
                    value = 0;

                    // Spin trace
                    value += reducedTwoParticleOperator(2*p, 2*q, 2*r, 2*s);
                    value += reducedTwoParticleOperator(2*p+1, 2*q+1, 2*r+1, 2*s+1);

                    value += reducedTwoParticleOperator(2*p+1, 2*q, 2*r, 2*s+1);
                    value += reducedTwoParticleOperator(2*p, 2*q+1, 2*r+1, 2*s);

                    rho2.insert( pair<int,cx_double>(mapTwoParticleStates(p,q,r,s), value) );
                }
            }
        }
    }

#ifdef DEBUG
//#if 1
//    cout << "n rho2 = " << rho2.size() << endl;
//    cout << "OrbitalEquation::computeTwoParticleReducedDensity()" << endl;
    // Taking the trace
    cx_double tra = 0;
    for(int i=0; i<nOrbitals; i++){
        for(int j=0; j<nOrbitals; j++){
            tra += findRho2(i, j, j, i);
        }
    }
    cout << "Trace(rho2) = " << real(tra) << " N(N-1) = " << nParticles*(nParticles-1) << endl;
#endif
}
//------------------------------------------------------------------------------
cx_double OrbitalEquation::reducedOneParticleOperator(const int i, const int j)
{
    bitset<BITS> newState;
    double phase;
    cx_double value = 0;

    for(int n=0; n<nSlaterDeterminants; n++){
        for(int m=0; m<nSlaterDeterminants; m++){
            phase = 1;
            newState = slaterDeterminants[m];

            //------------------------------------------------------------------
            // The following is an ugly way of writing
            //------------------------------------------------------------------
            //            removeParticle(j, newState);
            //            phase *= sign(j, newState);

            //            addParticle(i, newState);
            //            phase *= sign(i, newState);

            //            if(slaterDeterminants[n] == newState)
            //                value += phase*conj((*A)(n))*(*A)(m);
            //------------------------------------------------------------------

            removeParticle(j, newState);
            if(newState[BITS-1] != 1){
                phase *= sign(j, newState);

                addParticle(i, newState);

                if(newState[BITS-1] != 1){
                    phase *= sign(i, newState);

                    if (!newState[BITS - 1]) {
                        if(slaterDeterminants[n] == newState)
                            value += phase*conj((*A)(n))*(*A)(m);
                    }
                }
            }
            //------------------------------------------------------------------
        }
    }
    return value;
}
//------------------------------------------------------------------------------
cx_double OrbitalEquation::reducedTwoParticleOperator(const int p, const int q,
                                                      const int r, const int s)
{

    bitset<BITS> newState;
    double phase;
    cx_double value = 0;

    // Summing over all Slater determinants
    for(int n=0; n<nSlaterDeterminants; n++){
        for(int m=0; m<nSlaterDeterminants; m++){
            phase = 1;
            newState = slaterDeterminants[m];

            //------------------------------------------------------------------
            // The following is an ugly way of writing
            //------------------------------------------------------------------
            //            removeParticle(s, newState);
            //            phase *= sign(s, newState);

            //            removeParticle(r, newState);
            //            phase *= sign(r, newState);

            //            addParticle(q, newState);
            //            phase *= sign(q, newState);

            //            addParticle(p, newState);
            //            phase *= sign(p, newState);

            //            if(slaterDeterminants[n] == newState)
            //                value += phase*conj((*A)(n))*(*A)(m);
            //------------------------------------------------------------------
            removeParticle(s, newState);
            if(newState[BITS-1] != 1){
                phase *= sign(s, newState);

                removeParticle(r, newState);
                if(newState[BITS-1] != 1){
                    phase *= sign(r, newState);

                    addParticle(q, newState);
                    if(newState[BITS-1] != 1){
                        phase *= sign(q, newState);

                        addParticle(p, newState);
                        if(newState[BITS-1] != 1){
                            phase *= sign(p, newState);

                            if (!newState[BITS - 1]) {
                                if(slaterDeterminants[n] == newState)
                                    value += phase*conj((*A)(n))*(*A)(m);
                            }
                        }
                    }
                }
            }
            //------------------------------------------------------------------
        }
    }

    return value;
}
//------------------------------------------------------------------------------
#define WITH_SPIN 0
const cx_mat &OrbitalEquation::reCalculateRho1(const cx_vec &A)
{
    this->A = &A;
#if WITH_SPIN
    computeOneParticleReducedDensityWithSpin();
    return rho1;
#else
    computeOneParticleReducedDensity();
    return invRho; // This is not actually the inverse at this point
#endif
}
//------------------------------------------------------------------------------
vec OrbitalEquation::getSvdRho1(){
    cx_mat X;
    vec s;
    cx_mat Y;
#if WITH_SPIN
    svd_econ(X, s, Y, rho1);
#else
    svd_econ(X, s, Y, invRho);
#endif
    return s;
}
//------------------------------------------------------------------------------
double OrbitalEquation::getCorrelation()
{
    cx_mat rho1_sq;
#if WITH_SPIN
    rho1_sq = rho1*rho1/(pow(nParticles,2));
    cout << trace(rho1_sq) << endl;
    return 1.0/real(trace(rho1_sq));
#else
//    computeOneParticleReducedDensity();
    rho1_sq = invRho*invRho/(nParticles*nParticles);
#endif
    return 1.0/real(trace(rho1_sq));
}
//------------------------------------------------------------------------------
cx_double OrbitalEquation::findRho2(int p, int q, int r, int s)
{
    cx_double result = 0;
    int searchValue = mapTwoParticleStates(p,q,r,s);

    auto foundInteraction = rho2.find(searchValue);
    if(foundInteraction != rho2.end())
        result = foundInteraction->second;
    return result;
}
//------------------------------------------------------------------------------
void OrbitalEquation::computeOneParticleReducedDensityWithSpin()
{
    for(int i=0; i < 2*nOrbitals; i++){
        rho1(i,i) = reducedOneParticleOperator(i,i);
        for(int j=i+1; j < 2*nOrbitals; j++){
            rho1(i,j) = reducedOneParticleOperator(i,j);
            rho1(j,i) = conj( rho1(i,j) );
        }
    }

}
//------------------------------------------------------------------------------
