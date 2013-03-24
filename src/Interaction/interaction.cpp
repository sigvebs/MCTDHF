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
#if 0
    double threshold = 1e-10;

    cx_double V;
    for (int p = 0; p < nOrbitals; p++) {
        for (int q = 0; q < nOrbitals; q++) {
            for (int r = 0; r < nOrbitals; r++) {
                for (int s = 0; s < nOrbitals; s++) {
                    V = mfIntegrator->integrate(p,q,r,s, C);
#if 0
                    if(abs(real(V)) < threshold){
                        V = cx_double(0, imag(V));
                    }
                    if(abs(imag(V)) < threshold){
                        V = cx_double(real(V), 0);
                    }
                    //                    cout << p << q << r << s;
                    //                    cout << " V = " << V << endl;
                    //                    cout << "V = " << V << endl;
#endif


                    interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(p,q,r,s), V) );
                }
            }
        }
    }
#endif
#if 1
    cx_double V;
    for (int p = 0; p < nOrbitals; p++) {
        for (int q = p; q < nOrbitals; q++) {
            for (int r = 0; r < nOrbitals; r++) {
                for (int s = r; s < nOrbitals; s++) {
                    V = mfIntegrator->integrate(p,q,r,s, C);
                    interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(p,q,r,s), V) );
                    interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(q,p,s,r), V) );

                    V = mfIntegrator->integrate(p,q,s,r, C);
                    interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(p,q,s,r), V) );
                    interactionElements.insert( pair<int,cx_double>(mapTwoParticleStates(q,p,r,s), V) );
                }
            }
        }
    }
#endif
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
    delete mfIntegrator;
}
//------------------------------------------------------------------------------
void Interaction::printInteractionElements(){
    cout << "Interaction elements" << endl;
    for(auto it = interactionElements.begin(); it != interactionElements.end(); it++){
        cx_double a = it->second;
        cout << "id = " << it->first << " value = \t" << real(a) << "\t" << imag(a) << endl;
    }
}
//------------------------------------------------------------------------------
