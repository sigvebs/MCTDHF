#ifndef SLATERDETERMINANTS_H
#define SLATERDETERMINANTS_H


#include "src/includes/defines.h"

#include <cstdlib>
//#include <stdio.h>
//#include <iostream>
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
    const vector<bitset<BITS> > &getSlaterDeterminants() const;
    const vector<bool*> &getSlaterDeterminantsBool() const;
protected:
    vec odometer(const vec &, int, int);
    bitset<BITS> createBinaryState(vec state);
    bool checkEigenSpin(vec state);


    // Bool operations
    bool *createBinaryStateBool(vec state);

    Config *cfg;

    // System
    int nParticles;
    bool conservedEigenSpin;
    int conservedEigenSpinValue;

    vector<vec> sps;
    vector<bitset<BITS> > binStates;
    vector<bool*> binStatesBool;
};

#endif // SLATERDETERMINANTS_H
