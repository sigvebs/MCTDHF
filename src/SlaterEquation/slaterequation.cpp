#include "slaterequation.h"
//------------------------------------------------------------------------------
SlaterEquation::SlaterEquation(Config *cfg,
                               vector<bitset<BITS> > slaterDeterminants,
                               Interaction* interaction, SingleParticleOperator* oneParticleOperator):
    cfg(cfg),
    interaction(interaction),
    oneParticleOperator(oneParticleOperator)
{
    this->slaterDeterminants = slaterDeterminants;
    this->nSlaterDeterminants = slaterDeterminants.size();
    try {
        dim = cfg->lookup("system.dim");
        nOrbitals = cfg->lookup("spatialDiscretization.nSpatialOrbitals");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "SlaterEquation::Error reading from config object." << endl;
        exit(EXIT_FAILURE);
    }
    H = zeros<cx_mat>(nSlaterDeterminants, nSlaterDeterminants);
}
//------------------------------------------------------------------------------
cx_vec SlaterEquation::computeRightHandSide(const cx_vec &A)
{
    computeHamiltonianMatrix();

    return H*A;
}
//------------------------------------------------------------------------------
cx_vec SlaterEquation::computeRightHandSideComplexTime(const cx_vec &A)
{
    cout.precision(16);
    computeHamiltonianMatrix();

    cx_vec HA = H*A;
    cx_double E = cdot(A, HA)/cdot(A,A);

    return HA - E*A;
}
//------------------------------------------------------------------------------
void SlaterEquation::computeHamiltonianMatrix()
{
    cx_double phase;
    cx_double Vpqrs;
    h = &oneParticleOperator->getH();
    H.zeros();

    for(int m=0; m<nSlaterDeterminants; m++){
        for(int n=m; n<nSlaterDeterminants; n++){
            //------------------------------------------------------------------
            for(int p=0; p<nOrbitals; p++){
                for(int q=0; q<nOrbitals; q++){
                    //----------------------------------------------------------
                    // One body part

                    // The orbitals in the Slater determinants are mapped by
                    // 2*p == spin up
                    // 2*p+1 == spin down
                    // where p is the same spatial orbital

                    phase = 0;
                    phase += secondQuantizationOneBodyOperator(2*q, 2*p,
                                                                slaterDeterminants[n],
                                                                slaterDeterminants[m]);
                    phase += secondQuantizationOneBodyOperator(2*q+1, 2*p+1,
                                                                slaterDeterminants[n],
                                                                slaterDeterminants[m]);

                    H(m,n) += (*h)(p,q)*phase;

                    //----------------------------------------------------------

                    // Two-body interaction
                    for(int r=0; r<nOrbitals; r++){
                        for(int s=0; s<nOrbitals; s++){
                            Vpqrs = interaction->at(p, q, r, s);

                            if(Vpqrs != cx_double(0,0)){
                                phase = 0;


                                phase += secondQuantizationTwoBodyOperator(2*p, 2*q, 2*r, 2*s,
                                                                           slaterDeterminants[n],
                                                                           slaterDeterminants[m]);

                                phase += secondQuantizationTwoBodyOperator(2*p+1, 2*q+1, 2*r+1, 2*s+1,
                                                                           slaterDeterminants[n],
                                                                           slaterDeterminants[m]);

                                phase += secondQuantizationTwoBodyOperator(2*p, 2*q+1, 2*r, 2*s+1,
                                                                           slaterDeterminants[n],
                                                                           slaterDeterminants[m]);

                                phase += secondQuantizationTwoBodyOperator(2*p+1, 2*q, 2*r+1, 2*s,
                                                                           slaterDeterminants[n],
                                                                           slaterDeterminants[m]);

                                H(m, n) += 0.5*Vpqrs*phase;
                            }
                            //--------------------------------------------------
                        }
                    }

                    //----------------------------------------------------------
                }
            }
            //------------------------------------------------------------------
            if(m != n){
                H(n, m) = conj(H(m,n));
            }
        }
    }
#ifdef DEBUG
//#if 1
    cout << "SlaterEquation::computeHamiltonianMatrix()" << endl;
    //    cout << "Hamiltonian matrix:" << endl << H << endl;
    //    cx_vec eigval;
    //    cx_mat eigvec;
    //    eig_gen(eigval, eigvec, H);
    vec eigval = eig_sym(H);
    cout << "min(eigval) = " << min(eigval) << endl;
#endif
}
//------------------------------------------------------------------------------
cx_double SlaterEquation::secondQuantizationOneBodyOperator(const int p, const int q,
                                                           bitset<BITS> state1,
                                                           const bitset<BITS> &state2)
{
    cx_double phase = 1;

    removeParticle(q, state1);
    if(state1[BITS-1] == 1)
        return 0;
    phase *= sign(q, state1);

    addParticle(p, state1);
    if(state1[BITS-1] == 1)
        return 0;
    phase *= sign(p, state1);

    if (state2 !=  state1)
        phase = 0;

    return phase;
}

//------------------------------------------------------------------------------
cx_double SlaterEquation::secondQuantizationTwoBodyOperator(const int p, const int q,
                                                            const int r, const int s,
                                                           bitset<BITS> state1,
                                                           const bitset<BITS> &state2)
{
    cx_double phase = 1;
    removeParticle(r, state1);
    if(state1[BITS-1] == 1)
        return 0;
    phase *= sign(r, state1);

    removeParticle(s, state1);
    if(state1[BITS-1] == 1)
        return 0;
    phase *= sign(s, state1);

    addParticle(q, state1);
    if(state1[BITS-1] == 1)
        return 0;
    phase *= sign(q, state1);

    addParticle(p, state1);
    if(state1[BITS-1] == 1)
        return 0;
    phase *= sign(p, state1);

    if (state1 !=  state2)
        return 0;

    return phase;
}
//------------------------------------------------------------------------------
// This function is temporary for testing
double SlaterEquation::getEnergy(const cx_vec& A)
{
    computeHamiltonianMatrix();
    vec eigval = eig_sym(H);
    cout << "min(eigval) = " << min(eigval) << endl;

    return real(cdot(A, H*A)/cdot(A,A));
}
//------------------------------------------------------------------------------
SlaterEquation::~SlaterEquation()
{
}
//------------------------------------------------------------------------------
