#include "orbitalequation.h"

//------------------------------------------------------------------------------
OrbitalEquation::OrbitalEquation(Config *cfg,
                                 vector<vec> orbitals,
                                 vector<bitset<BITS> > slaterDeterminants,
                                 Interaction *V,
                                 SingleParticleOperator *h):
    cfg(cfg),
    orbitals(orbitals),
    slaterDeterminants(slaterDeterminants),
    V(V),
    h(h)
{
    double L;
    double dx;
    try {
        nParticles = cfg->lookup("system.nParticles");
        dim = cfg->lookup("system.dim");
        dx = cfg->lookup("spatialDiscretization.gridSpacing");
        L = cfg->lookup("spatialDiscretization.latticeRange");
        nGrid = cfg->lookup("spatialDiscretization.nGrid");
        nSpatialOrbitals = cfg->lookup("spatialDiscretization.nSpatialOrbitals");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "OrbitalEquation::Error reading from config object." << endl;
        exit(EXIT_FAILURE);
    }

    nOrbitals = orbitals.size()/2;
    nSlaterDeterminants = slaterDeterminants.size();
    rho = cx_mat(nSpatialOrbitals, nSpatialOrbitals);

    U = cx_mat(nGrid, nSpatialOrbitals);
    P = cx_mat(nGrid, nGrid);
    hC = cx_mat(nGrid, nSpatialOrbitals);
    I = eye<cx_mat>(nGrid, nGrid);
}
//------------------------------------------------------------------------------
cx_mat OrbitalEquation::computeRightHandSide(cx_mat C, cx_vec A)
{
    this->C = C;
    this->A = A;

    computeProjector();
    computeOneParticleReducedDensity();
    computeTwoParticleReducedDensity();
//    computeUMatrix();
    // Computing the right hand side of the equation
    hC = h->getHspatial();

    cx_mat rightHandSide = (I - P)*hC;
    //    cout << rightHandSide << endl;

//    rightHandSide = zeros<cx_mat>(nGrid, nOrbitals);
    return rightHandSide;
}
//------------------------------------------------------------------------------
void OrbitalEquation::computeUMatrix()
{
    cx_double V_C;
    int searchValue;
    double intElement;
    cx_double invRho;
    cx_double Utmp;

    for(int j=0; j<nSpatialOrbitals; j++){

        for(int x=0; x<nGrid; x++){
            U(x,j) = 0;
            Utmp = 0;
            for(int i=0; i<nSpatialOrbitals; i++){
                invRho = rho(j,i);
                for(int q=0; q<nSpatialOrbitals; q++){
                    for(int r=0; r<nSpatialOrbitals; r++){
                        for(int s=0; s<nSpatialOrbitals; s++){

                            searchValue = mapTwoParticleStates(i,q,s,r);
                            auto foundRho2 = rho2.find(searchValue);
                            if(foundRho2 != rho2.end()){
                                Utmp *= foundRho2->second;
                            }

                        }
                    }
                }
            }
        }
    }
    //    for(int i=0; i<nOrbitals; i++){
    //        for(int beta=0; beta<nConstStates; beta++){

    //            U(i, beta) = 0;

    //            for(int q=0; q<nOrbitals; q++){
    //                for(int r=0; r<nOrbitals; r++){
    //                    for(int s=0; s<nOrbitals; s++){

    //                        searchValue = mapTwoParticleStates(i,q,s,r);
    //                        auto foundRho2 = rho2.find(searchValue);

    //                        if(foundRho2 != rho2.end()){
    //                            V_C = 0;
    //
    //                                }
    //                            }
    //                            U(i, beta) += foundRho2->second*V_C;
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //    }

#ifdef DEBUG
    cout << "OrbitalEquation::computeUMatrix()" << endl;
    //        cout << "U = " << endl << U << endl;
    //        exit(0);
#endif
}
//------------------------------------------------------------------------------
void OrbitalEquation::computeProjector()
{
    P = zeros<cx_mat>(nGrid, nGrid);

    // Slightly changing the projector to ensure orthonormality is conserved
    // Problems arise due to numerical inaccuracies.
#if 0
    cx_mat O = zeros<cx_mat>(nSpatialOrbitals,nSpatialOrbitals);

    for(int i=0; i<nSpatialOrbitals; i++){
        for(int l=0; l<nSpatialOrbitals; l++){
            O(i,l) = cdot(C.col(i), C.col(l));
        }
    }
    cout << O << endl;
    O = inv(O);

    for(int i=0; i<nSpatialOrbitals; i++){
        for(int l=0; l<nSpatialOrbitals; l++){
            O = cdot(C.col(i), C.col(l));
            P += C.col(i)*C.col(l).t()*O(i,l);
        }
    }
#else
    // The ecaxt mathematical defintion of the projector.
    for(int i=0; i<nSpatialOrbitals; i++){
        P += C.col(i)*C.col(i).t();
    }
#endif
#ifdef DEBUG
    cout << "OrbitalEquation::computeProjector()" << endl;
#endif
}
//------------------------------------------------------------------------------
void OrbitalEquation::computeOneParticleReducedDensity()
{
    bitset<BITS> newState;
    double phase;
    // All one particle combinantions
    for(int i=0; i< nSpatialOrbitals; i++){
        for(int j=i; j<nSpatialOrbitals; j++){
            rho(i,j) = 0;

            // Spin up
            for(int n=0; n<nSlaterDeterminants; n++){
                for(int m=0; m<nSlaterDeterminants; m++){
                    phase = 1;
                    newState = slaterDeterminants[m];

                    removeParticle(2*j, newState);
                    phase *= sign(2*j, newState);

                    addParticle(2*i, newState);
                    phase *= sign(2*i, newState);

                    if (!newState[BITS - 1]) {
                        if(slaterDeterminants[n] == newState)
                            rho(i,j) += phase*conj(A(n))*A(m);
                    }
                }
            }

            // Spin down
            for(int n=0; n<nSlaterDeterminants; n++){
                for(int m=0; m<nSlaterDeterminants; m++){
                    phase = 1;
                    newState = slaterDeterminants[m];

                    removeParticle(2*j + 1, newState);
                    phase *= sign(2*j + 1, newState);

                    addParticle(2*i +1, newState);
                    phase *= sign(2*i + 1, newState);

                    if (!newState[BITS - 1]) {
                        if(slaterDeterminants[n] == newState)
                            rho(i,j) += phase*conj(A(n))*A(m);
                    }
                }
            }
            rho(j,i) = conj(rho(i,j));
        }
    }

    rho = inv(rho);
#ifdef DEBUG
    cout << "OrbitalEquation::computeOneParticleReducedDensity()" << endl;
    cout << "Trace(rho) = " << real(trace(rho)) << " nParticles = " << nParticles << endl;
#endif
}
//------------------------------------------------------------------------------
void OrbitalEquation::computeTwoParticleReducedDensity()
{
    rho2.clear();
    cx_double value;

    // All two particle combinantions
    for(int p=0; p<nSpatialOrbitals; p++){
        for(int q=0; q<nSpatialOrbitals; q++){
            for(int r=0; r<nSpatialOrbitals; r++){
                for(int s=0; s<nSpatialOrbitals; s++){
                    value = 0;

                    // Spin trace
                    value += reducedTwoParticleOperator(2*p, 2*q, 2*r, 2*s);
                    value += reducedTwoParticleOperator(2*p+1, 2*q+1, 2*r+1, 2*s+1);

                    value += reducedTwoParticleOperator(2*p, 2*q+1, 2*r, 2*s+1);
                    value += reducedTwoParticleOperator(2*p+1, 2*q, 2*r+1, 2*s);

                    if(real(conj(value)*value) > 0){
                        rho2.insert( pair<int,cx_double>(mapTwoParticleStates(p,q,r,s), value) );
                        //                        rho2.insert( pair<int,cx_double>(mapTwoParticleStates(q,p,s,r), value) );
                        //                        rho2.insert( pair<int,cx_double>(mapTwoParticleStates(p,q,s,r), -value) );
                        //                        rho2.insert( pair<int,cx_double>(mapTwoParticleStates(q,p,r,s), -value) );
                    }
                }
            }
        }
    }

#ifdef DEBUG
    cout << "OrbitalEquation::computeTwoParticleReducedDensity()" << endl;
    // Taking the trace
    cx_double trace = 0;
    for(int i=0; i<nSpatialOrbitals; i++){
        for(int j=0; j<nSpatialOrbitals; j++){
            int searchValue = mapTwoParticleStates(i,j,i,j);
            auto found = rho2.find(searchValue);
            if(found != rho2.end()){
                trace += found->second;
            }
        }
    }
    cout << "Trace(rho2) = " << trace << " nParticles = " << nParticles << endl;
#endif
}
//------------------------------------------------------------------------------
cx_double OrbitalEquation::reducedTwoParticleOperator(int p, int q, int s, int r)
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
                if(slaterDeterminants[n] == newState){
                    value += phase*conj(A(n))*A(m);
                }
            }
        }
    }

    return value;
}
//------------------------------------------------------------------------------
void OrbitalEquation::setInititalState(cx_mat C)
{
    //    cx_mat X = randu<cx_mat>(nConstStates, nOrbitals);
    //    cx_mat R;
    //    qr(C,R,X);

    //    C = eye<cx_mat>(nGrid, nOrbitals);
    this->C = C;
    //        cout << "C = " << C << endl;

#ifdef DEBUG
    cout << "const cx_mat &OrbitalEquation::setInititalState()" << endl;
#endif
}
//------------------------------------------------------------------------------
void OrbitalEquation::computeOneBodyMatrix()
{
    hC = h->getH()*C;
#ifdef DEBUG
    cout << "OrbitalEquation::computeOneBodyMatrix()" << endl;
    //    cout << hC << endl;
    //    cout << C << endl;
#endif
}
//------------------------------------------------------------------------------
void OrbitalEquation::setSlaterCoefficients(const cx_vec &A)
{
    this->A = A;
}
//------------------------------------------------------------------------------
const cx_mat &OrbitalEquation::getCoefficientMatrix() const
{
    return C;
}
//------------------------------------------------------------------------------
