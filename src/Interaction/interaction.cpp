#include "interaction.h"
//------------------------------------------------------------------------------
Interaction::Interaction(Config *cfg, MeanFieldIntegrator *mfIntegrator):
    cfg(cfg),
    mfIntegrator(mfIntegrator)
{
    try{
        nOrbitals = cfg->lookup("spatialDiscretization.nSpatialOrbitals");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "Interaction::Interaction(Config *cfg)"
             << "::Error reading from config object." << endl;
    }
}

//------------------------------------------------------------------------------
void Interaction::updatePositionBasisElements(){
    mfIntegrator->initialize();
}
//------------------------------------------------------------------------------
void Interaction::computeNewElements(const cx_mat &C)
{
    interactionElements.clear();

    mfIntegrator->computeMeanField(C);
    computeInteractionelements(C);
}
//------------------------------------------------------------------------------
void Interaction::computeInteractionelements(const cx_mat &C)
{
    double tolerance = 1e-13;

    cx_double V;
    for (int p = 0; p < nOrbitals; p++) {
        for (int q = 0; q < nOrbitals; q++) {
            for (int r = 0; r < nOrbitals; r++) {
                for (int s = 0; s < nOrbitals; s++) {
                    V = mfIntegrator->integrate(p,q,r,s, C);

//                    if(abs(real(V)) < tolerance){
//                        V = cx_double(0, imag(V));
//                    }

//                    V = cx_double(real(V), 0);
//                    if(abs(imag(V)) < tolerance){
//                        V = cx_double(real(V), 0);
//                    }

//                    if (abs(V) > tolerance)
                        interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(p,q,r,s), V) );


                }
            }
        }
    }
//    cx_double V;
//    for (int p = 0; p < nOrbitals; p++) {
//        for (int q = p; q < nOrbitals; q++) {
//            for (int r = 0; r < nOrbitals; r++) {
//                for (int s = r; s < nOrbitals; s++) {
//                    // Symmetric
//                    V = mfIntegrator->integrate(p,q,r,s, C);
////                    V = cx_double(real(V), 0);

////                    if (abs(V) > tolerance) {
//                        interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(p,q,r,s), V) );
//                        interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(q,p,s,r), V) );
////                        interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(r,s,p,q), conj(V)) );
////                        interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(s,r,q,p), conj(V)) );
////                    }

//                    // Anti-Symmetric
//                    V = mfIntegrator->integrate(p,q,s,r, C);
////                    V = cx_double(real(V), 0);

////                    if (abs(V) > tolerance) {
//                        interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(p,q,s,r), V) );
//                        interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(q,p,r,s), V) );
////                        interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(s,r,p,q), conj(V)) );
////                        interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(r,s,q,p), conj(V)) );
////                    }
//                }
//            }
//        }
//    }

#ifdef DEBUG
//#if 1
//    cout << "Number of interaction elements = " << interactionElements.size() << endl;
//    for(auto it = interactionElements.begin(); it != interactionElements.end(); it++){
//        cout << "id = " << it->first << " value = " << it->second << endl;
//    }
//    cout << "nInteraction = " << interactionElements.size() << endl;
//        exit(1);
#endif
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
void Interaction::addPotential(InteractionPotential *interactionPot)
{
    mfIntegrator->addPotential(interactionPot);
}
//------------------------------------------------------------------------------
Interaction::~Interaction()
{

}
//------------------------------------------------------------------------------
