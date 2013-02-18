#include "interaction.h"
//------------------------------------------------------------------------------
Interaction::Interaction(Config *cfg, vector<vec> orbitals):
    cfg(cfg),
    orbitals(orbitals)
{

    double L;
    try{
        L = cfg->lookup("spatialDiscretization.latticeRange");
        dx = cfg->lookup("spatialDiscretization.gridSpacing");
        aa = cfg->lookup("potential.a");
        nGrid = cfg->lookup("spatialDiscretization.nGrid");
        nSpatialOrbitals = cfg->lookup("spatialDiscretization.nSpatialOrbitals");
        aa *=aa;
    } catch (const SettingNotFoundException &nfex) {
        cerr << "BasisHarmonicOscillator::BasisHarmonicOscillator(Config *cfg)"
             << "::Error reading from config object." << endl;
    }

    nOrbitals = orbitals.size();
    V2 = field<cx_vec>(nSpatialOrbitals, nSpatialOrbitals);
    x = linspace(-L, L, nGrid);
}
//------------------------------------------------------------------------------
void Interaction::computeNewElements(const cx_mat &C)
{
    this->C = C;
    computeMeanField();
    computeInteractionelements();
}
//------------------------------------------------------------------------------
void Interaction::computeMeanField()
{
    for (int q = 0; q < nSpatialOrbitals; q++) {
        for (int r = q; r < nSpatialOrbitals; r++) {
            V2(q,r) = integrate(q,r);
            V2(r,q) = conj(V2(q,r));
        }
    }
}
//------------------------------------------------------------------------------
void Interaction::computeInteractionelements()
{
    double tolerance = 1e-6;

    cx_double V;
    for (int p = 0; p < nOrbitals; p++) {
        for (int q = p + 1; q < nOrbitals; q++) {
            for (int r = 0; r < nOrbitals; r++) {
                for (int s = r + 1; s < nOrbitals; s++) {

                    // Symmetric spin
                    V = 0;
                    if (orbitals[p][0] == orbitals[r][0] && orbitals[q][0] == orbitals[s][0])
                        V = integrate(p,q,r,s);

                    if (abs(real(V)) > tolerance) {
//                        cout << p << q << r << s << " V = " << real(V) << endl;
                        interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(p,q,r,s), V) );
                        interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(q,p,s,r), conj(V)) );
                    }

                    // Anti-Symmetric spin
                    V = 0;
                    if (orbitals[p][0] == orbitals[s][0] && orbitals[q][0] == orbitals[r][0])
                        V = integrate(p,q,s,r);

                    if (abs(real(V)) > tolerance) {
//                        cout << p << q << s << r << " V = " << real(V) << endl;
                        interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(p,q,s,r), V) );
                        interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(q,p,r,s), conj(V)) );
                    }
                }
            }
        }
    }

//     TODO: REMOVE ME!!!!
    interactionElements.clear();

#ifdef DEBUG
    cout << "void Basis::computeInteractionelements()" << endl;
//    cout << "Number of interaction elements = " << interactionElements.size() << endl;
//        for(auto it = interactionElements.begin(); it != interactionElements.end(); it++){
//            cout << "id = " << it->first << " value = " << it->second << endl;
//        }
    //        int searchValue = mapTwoParticleStates(0,2,0,2);
    //        auto found = interactionElements.find(searchValue);

    //        if(found != interactionElements.end() )
    //            cout << "found = " << found->first << " " << found->second << endl;
    //        else
    //            cout << "not found" << endl;
#endif
}
//------------------------------------------------------------------------------
cx_vec Interaction::integrate(int q, int r)
{
    cx_vec V = zeros<cx_vec>(nGrid);
    cx_double integral ;

    for(int i=0; i<nGrid; i++){
        integral = 0;
        for(int j=0; j<nGrid; j++){
            integral += conj(C(j,q))/sqrt(pow(x[j] - x[i],2) + aa)*C(j,r);
        }
        V(i) = integral;
    }

    return V;
}
//------------------------------------------------------------------------------
cx_double Interaction::integrate(int p, int q, int r, int s)
{
    cx_double integral = 0;

    for(int i=0; i<nGrid; i++){
//        integral += conj(C(i,p))*V2(q,s)(i)*C(i,r);
        integral += conj(C(i,p/2))*V2(q/2,s/2)(i)*C(i,r/2);
    }

    return integral;
}
//------------------------------------------------------------------------------
cx_double Interaction::at(const int p, const int q, const int r, const int s)
{
    cx_double result = 0;
    int searchValue = mapTwoParticleStates(p,q,r,s);

    auto foundInteraction = interactionElements.find(searchValue);
    if(foundInteraction != interactionElements.end())
        result = foundInteraction->second;
    return result;
}
//------------------------------------------------------------------------------
