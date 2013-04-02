#ifndef SLATERDETERMINANTS_H
#define SLATERDETERMINANTS_H

// Local includes
#include "src/includes/defines.h"

// Library includes
#include <cstdlib>
#include <armadillo>
#include <vector>
#include <bitset>
#include <libconfig.h++>

using namespace libconfig;
using namespace std;
using namespace arma;

class SlaterDeterminants
{
public:
    SlaterDeterminants(Config *cfg, vector<vec> sps);
    void createSlaterDeterminants();
    void createInitialState();
    const vector<bitset<BITS> > &getSlaterDeterminants() const;
    void saveSlaterDeterminantsToDisk();
    cx_vec getCoefficients();
    void load();
protected:
    vector<vec> sps;
    vector<bitset<BITS> > binStates;
    cx_vec A;

    vec odometer(const vec &, int, int);
    bitset<BITS> createBinaryState(vec state);
    bool checkEigenSpin(vec state);

    // System
    int nParticles;
    bool conservedEigenSpin;
    int conservedEigenSpinValue;
    string filePath;

    Config *cfg;
};

#endif // SLATERDETERMINANTS_H
