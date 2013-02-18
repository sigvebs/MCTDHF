#include "slaterequation.h"
//------------------------------------------------------------------------------
SlaterEquation::SlaterEquation(Config *cfg,
                               vector<vec> orbitals,
                               vector<bitset<BITS> > slaterDeterminants,
                               Interaction* interaction, SingleParticleOperator* oneParticleOperator):
    cfg(cfg),
    orbitals(orbitals),
    nOrbitals(orbitals.size()),
    slaterDeterminants(slaterDeterminants),
    nSlaterDeterminants(slaterDeterminants.size()),
    interaction(interaction),
    oneParticleOperator(oneParticleOperator)
{
    try {
        dim = cfg->lookup("system.dim");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "SlaterEquation::Error reading from config object." << endl;
        exit(EXIT_FAILURE);
    }
}
//------------------------------------------------------------------------------
cx_vec SlaterEquation::computeRightHandSide(const cx_vec &A)
{
    this->A = A;
    h = oneParticleOperator->getH();

    computeHamiltonianMatrix();
    cx_vec HA = H*A;

    return HA;
}
//------------------------------------------------------------------------------
cx_vec SlaterEquation::computeRightHandSideComplexTime(const cx_vec &A)
{
    this->A = A;
    h = oneParticleOperator->getH();

    computeHamiltonianMatrix();
    cx_vec HA = H*A;
    cx_double E = cdot(A, HA)/cdot(A,A);

    return HA - E*A;
}
//------------------------------------------------------------------------------
void SlaterEquation::computeHamiltonianMatrix()
{
    cx_double phase;
    bitset<BITS> newState;

    for(int m=0; m<nSlaterDeterminants; m++){
        for(int n=m; n<nSlaterDeterminants; n++){
            H(m,n) = 0;
            //------------------------------------------------------------------
            for(int p=0; p<nOrbitals; p++){
                for(int q=0; q<nOrbitals; q++){
                    //----------------------------------------------------------
                    // One body
                    phase = 1;
                    newState = slaterDeterminants[n];

                    removeParticle(q, newState);
                    phase *= sign(q, newState);

                    addParticle(p, newState);
                    phase *= sign(p, newState);

                    if (newState ==  slaterDeterminants[m]){
                        H(m,n) += h(p,q)*phase;
                    }

                    //----------------------------------------------------------
                    // Two-body interaction
                    for(int r=0; r<nOrbitals; r++){
                        for(int s=0; s<nOrbitals; s++){

                            phase = 1;
                            newState = slaterDeterminants[n];

                            removeParticle(r, newState);
                            phase *= sign(r, newState);

                            removeParticle(s, newState);
                            phase *= sign(s, newState);

                            addParticle(q, newState);
                            phase *= sign(q, newState);

                            addParticle(p, newState);
                            phase *= sign(p, newState);

                            if (newState ==  slaterDeterminants[m]){
                                H(m, n) += 0.5*phase*interaction->at(p, q, r, s);
                            }
                            //--------------------------------------------------
                        }
                    }
                    //----------------------------------------------------------
                }
            }
            //------------------------------------------------------------------
            H(n, m) = conj(H(m,n));
        }
    }
#ifdef DEBUG
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
void SlaterEquation::setInitalState()
{
    A = randu<cx_vec>(nSlaterDeterminants);
    A = A/sqrt(cdot(A, A));
    H = zeros<cx_mat>(nSlaterDeterminants, nSlaterDeterminants);
}
//------------------------------------------------------------------------------
const cx_mat &SlaterEquation::getCoefficientVector() const
{
    return A;
}
//------------------------------------------------------------------------------
// This function is temporary for testing
double SlaterEquation::getEnergy()
{
    h = oneParticleOperator->getH();
    computeHamiltonianMatrix();
    vec eigval = eig_sym(H);

    return real(cdot(A, H*A));
//    return min(eigval);
}
//------------------------------------------------------------------------------
SlaterEquation::~SlaterEquation()
{
}
//------------------------------------------------------------------------------
