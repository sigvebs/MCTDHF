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
        dim = cfg->lookup("system.dim");
        dx = cfg->lookup("spatialDiscretization.gridSpacing");
        nGrid = cfg->lookup("spatialDiscretization.nGrid");
        nOrbitals = cfg->lookup("spatialDiscretization.nSpatialOrbitals");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "OrbitalEquation::Error reading from config object." << endl;
        exit(EXIT_FAILURE);
    }

    nSlaterDeterminants = slaterDeterminants.size();
    invRho = cx_mat(nOrbitals, nOrbitals);

    U = zeros<cx_mat>(nGrid, nOrbitals);
    Q = cx_mat(nGrid, nGrid);
}
//------------------------------------------------------------------------------
const cx_mat &OrbitalEquation::computeRightHandSide(const cx_mat &C, const cx_vec &A)
{
    this->A = &A;
    hC = &(h->getHspatial());

    // Clearing values
    rho2.clear();

    computeProjector(C);
    computeOneParticleReducedDensity();
    computeTwoParticleReducedDensity();
    computeUMatrix(C);

    // Computing the right hand side of the equation
    rightHandSide = Q*U;
    return rightHandSide;
}
//------------------------------------------------------------------------------
void OrbitalEquation::computeUMatrix(const cx_mat &C)
{
    cx_vec Ui(nGrid);
    cx_vec inner(nGrid);
#if 1
    for(int j=0; j<nOrbitals; j++){
        U.col(j) = hC->col(j);

        for(int i=0; i<nOrbitals; i++){
//            Ui = zeros<cx_vec>(nGrid);
            for(int d=0; d<nGrid; d++)
                Ui(d) = 0;

            for(int q=0; q<nOrbitals; q++){
                for(int r=0; r<nOrbitals; r++){
//                    inner = zeros<cx_vec>(nGrid);
                    for(int d=0; d<nGrid; d++)
                        inner(d) = 0;

                    for(int s=0; s<nOrbitals; s++){
                        auto rho_iqrs = rho2.find(mapTwoParticleStates(i,q,r,s));
                        inner += rho_iqrs->second* C.col(s);
                    }
                    Ui += V->meanField(q,r)%inner;
//                    Ui += diagmat(V->meanField(q,r))*inner;
                }
            }
            U.col(j) += invRho(j,i)*Ui;
        }
    }
#else
    for(int j=0; j<nOrbitals; j++){
        U.col(j) = hC->col(j);
//        U.col(j) = zeros<cx_vec>(nGrid);
        Ui = zeros<cx_vec>(nGrid);
        for(int i=0; i<nOrbitals; i++){

            for(int q=0; q<nOrbitals; q++){
                for(int r=0; r<nOrbitals; r++){
                    inner = zeros<cx_vec>(nGrid);
                    for(int s=0; s<nOrbitals; s++){
                        auto rho_iqrs = rho2.find(mapTwoParticleStates(i,q,r,s));

//                        cout << " q = " << q << " r = " << " s = " << s
//                             << " rho2 = " << rho_iqrs->second << endl;
                        Ui += invRho(j,i)*rho_iqrs->second* diagmat(V->meanField(q,r))* C.col(s);
                    }
                }
            }
        }
        U.col(j) += Ui;
    }
#endif
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
#if 0
    // NEED TO BE UPDATED - not working
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
    P = -P; // TMP solution

    // OLD VERSION
    //    Q = zeros<cx_mat>(nGrid, nGrid);
    // The exact mathematical defintion of the projector.
    //    cx_mat Qn = zeros<cx_mat>(nGrid, nGrid);
    //    for(int i=0; i<nOrbitals; i++){
    //
    //    cout << max(max(Qn - Q));
    //    exit(1); Q -= C.col(i)*C.col(i).t();
    //    }
#else
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
    // First rho is calculated

    // All possible spatial orbitals
    for(int i=0; i< nOrbitals; i++){
        for(int j=i; j<nOrbitals; j++){
            invRho(i,j)  = reducedOneParticleOperator(2*i,2*j);
            invRho(i,j)  += reducedOneParticleOperator(2*i+1,2*j+1);
            invRho(j,i) = conj(invRho(i,j));
        }
    }

#ifdef DEBUG
    cout << "OrbitalEquation::computeOneParticleReducedDensity()" << endl;
    cout << "Trace(rho) = " << real(trace(rho)) << " nParticles = " << nParticles << endl;
//    exit(0);
#endif
    invRho = inv(invRho);
}
//------------------------------------------------------------------------------
void OrbitalEquation::computeTwoParticleReducedDensity()
{
    cx_double value;

    // All possible two spatial orbital combinantions
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

//                    cout << p << q << r << s << " V = " << value << endl;

                    rho2.insert( pair<int,cx_double>(mapTwoParticleStates(p,q,r,s), value) );
                    //                        rho2.insert( pair<int,cx_double>(mapTwoParticleStates(q,p,s,r), value) );
                    //                        rho2.insert( pair<int,cx_double>(mapTwoParticleStates(p,q,s,r), -value) );
                    //                        rho2.insert( pair<int,cx_double>(mapTwoParticleStates(q,p,r,s), -value) );
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
            auto found = rho2.find( mapTwoParticleStates(i,j,j,i) );
            if(found != rho2.end()){
                tra += found->second;
            }
        }
    }
    cout << "Trace(rho2) = " << tra << " N(N-1) = " << nParticles*(nParticles-1) << endl;
//    cout << rho1 << endl;
//    cout << rho << endl;
//    exit(1);
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

            removeParticle(j, newState);
            phase *= sign(j, newState);

            addParticle(i, newState);
            phase *= sign(i, newState);

            if (!newState[BITS - 1]) {
                if(slaterDeterminants[n] == newState)
                    value += phase*conj((*A)(n))*(*A)(m);
            }
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

            removeParticle(s, newState);
            phase *= sign(s, newState);

            removeParticle(r, newState);
            phase *= sign(r, newState);

            addParticle(q, newState);
            phase *= sign(q, newState);

            addParticle(p, newState);
            phase *= sign(p, newState);

            if (!newState[BITS - 1]) {
                if(slaterDeterminants[n] == newState)
                    value += phase*conj((*A)(n))*(*A)(m);
            }
        }
    }

    return value;
}
//------------------------------------------------------------------------------
