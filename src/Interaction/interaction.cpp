#include "interaction.h"
//------------------------------------------------------------------------------
Interaction::Interaction(Config *cfg):
    cfg(cfg)
{
    double L;
    try{
        L = cfg->lookup("spatialDiscretization.latticeRange");
        dx = cfg->lookup("spatialDiscretization.gridSpacing");
        aa = cfg->lookup("potential.a");
        nGrid = cfg->lookup("spatialDiscretization.nGrid");
        nOrbitals = cfg->lookup("spatialDiscretization.nSpatialOrbitals");
        aa *=aa;
    } catch (const SettingNotFoundException &nfex) {
        cerr << "BasisHarmonicOscillator::BasisHarmonicOscillator(Config *cfg)"
             << "::Error reading from config object." << endl;
    }

    V2 = field<cx_vec>(nOrbitals, nOrbitals);
    x = linspace(-L, L, nGrid);
}
//------------------------------------------------------------------------------
void Interaction::computeNewElements(const cx_mat &C)
{
    this->C = C;
    interactionElements.clear();
    V2 = field<cx_vec>(nOrbitals, nOrbitals);

    computeMeanField();
    computeInteractionelements();
}
//------------------------------------------------------------------------------
void Interaction::computeMeanField()
{
    for (int q = 0; q < nOrbitals; q++) {
        for (int r = q; r < nOrbitals; r++) {
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
        for (int q = p; q < nOrbitals; q++) {
            for (int r = 0; r < nOrbitals; r++) {
                for (int s = r; s < nOrbitals; s++) {
                    // Symmetric
                    V = integrate(p,q,r,s);

                    if (abs(real(V)) > tolerance) {
                        interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(p,q,r,s), V) );
                        interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(q,p,s,r), conj(V)) );
                    }

                    // Anti-Symmetric
                    V = integrate(p,q,s,r);

                    if (abs(real(V)) > tolerance) {
                        interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(p,q,s,r), V) );
                        interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(q,p,r,s), conj(V)) );
                    }
                }
            }
        }
    }

#ifdef DEBUG
//    cout << "Number of interaction elements = " << interactionElements.size() << endl;
//    cout << "nInteraction = " << interactionElements.size() << endl;
//        for(auto it = interactionElements.begin(); it != interactionElements.end(); it++){
//            cout << "id = " << it->first << " value = " << it->second << endl;
//        }
    //        int searchValue = mapTwoParticleStates(0,2,0,2);
    //        auto found = interactionElements.find(searchValue);

    //        if(found != interactionElements.end() )
    //            cout << "found = " << found->first << " " << found->second << endl;
    //        else
    //            cout << "not found" << endl;
//    exit(1);
#endif
}
//------------------------------------------------------------------------------
cx_vec Interaction::integrate(const int q, const int r)
{
    cx_vec V = zeros<cx_vec>(nGrid);
    cx_double integral ;

    // Calulating the mean field V^qr
    for(int i=0; i<nGrid; i++){

        // Integrations using the trapezodial rule.
        integral = 0;
        for(int j=1; j<nGrid-1; j++){
            integral += conj(C(j,q))/sqrt(pow(x(j) - x(i),2) + aa)*C(j,r);
        }
        integral *=2;

        // Enpoints
        integral += conj(C(0,q))/sqrt(pow(x(0) - x(i),2) + aa)*C(0,r)
                + conj(C(nGrid-1,q))/sqrt(pow(x(nGrid-1) - x(i),2) + aa)*C(nGrid-1,r);

        V(i) = 0.5*integral;
    }

    return V;
}
//------------------------------------------------------------------------------
cx_double Interaction::integrate(const int p, const int q, const int r, const int s)
{
    // Integrations using the trapezodial rule.
    cx_double integral = 0;
    for(int i=1; i<nGrid-1; i++){
        integral += conj(C(i,p))*V2(q,s)(i)*C(i,r);
    }
    integral *= 2;

    // Enpoints
    integral += conj(C(0,p))*V2(q,s)(0)*C(0,r) + conj(C(nGrid-1,p))*V2(q,s)(nGrid-1)*C(nGrid-1,r);

    return 0.5*integral;
}
//------------------------------------------------------------------------------
const cx_double Interaction::at(const int p, const int q, const int r, const int s)
{
    cx_double result = 0;
    int searchValue = mapTwoParticleStates(p,q,r,s);

    auto foundInteraction = interactionElements.find(searchValue);
    if(foundInteraction != interactionElements.end())
        result = foundInteraction->second;
    return result;
}
//------------------------------------------------------------------------------
const cx_vec &Interaction::meanField(const int p, const int q)
{
    return V2(p,q);
}
//------------------------------------------------------------------------------
