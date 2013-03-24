#include "basis.h"
//------------------------------------------------------------------------------
Basis::Basis(Config *cfg):
    cfg(cfg)
{
    string filePath;
    double L;
    bool periodicBoundaries;
    try{
        dim = cfg->lookup("system.dim");
        coordinateType = cfg->lookup("system.coordinateType");
        cfg->lookupValue("systemSettings.filePath", filePath);
        L = cfg->lookup("spatialDiscretization.latticeRange");
        nGrid = cfg->lookup("spatialDiscretization.nGrid");
        nBasis = cfg->lookup("system.shells");
        periodicBoundaries = cfg->lookup("spatialDiscretization.periodicBoundaries");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "Basis(Config *cfg)::Error reading from config object." << endl;
        exit(EXIT_FAILURE);
    }

    nSpatialOrbitals = states.size()/2;

    // Adding the number of gridpoints to the config file
    Setting &root = cfg->getRoot();
    Setting &tmp = root["spatialDiscretization"];

    tmp.add("gridSpacing", Setting::TypeFloat) = dx;

    if(periodicBoundaries){
        dx = 2.0*L/(double)(nGrid);
        cout << "hei" << endl;
        x = linspace<vec>(-L, L-dx,nGrid);
    }else{
        cout << "nei" << endl;
        x = linspace<vec>(-L, L, nGrid);
        dx = x(1) - x(0);
    }
    // Saving the grid basis
    filnameAxis = filePath + "x.mat";
    x.save(filnameAxis, arma_ascii);

#ifdef DEBUG
    cout << "Basis::Basis(Config *cfg)" << endl
         << "coordinateType \t= " << coordinateType << endl
         << "dim \t\t= " << dim << endl;
#endif
}
//------------------------------------------------------------------------------
void Basis::createBasis()
{
    switch(coordinateType){
    case CARTESIAN:
        createCartesianBasis();
        break;
    case POLAR:
        createPolarBasis();
        break;
    }

    // Adding the number of spatial orbitals to the config object
    nSpatialOrbitals = states.size()/2;
    Setting &root = cfg->getRoot();
    Setting &tmp = root["spatialDiscretization"];
    tmp.add("nSpatialOrbitals", Setting::TypeInt) = nSpatialOrbitals;
}
//------------------------------------------------------------------------------
void Basis::createCartesianBasis()
{
    // Generating all single particle states and energies
    vec state = zeros(dim + 1);

    switch (dim) {
    case 1:
        for (int n = 0; n <= nBasis; n++) {
            for (int spin = 1; spin >= -1; spin -= 2) {
                state(0) = spin;
                state(1) = n;
                states.push_back(state);
            }
        }
        break;
    case 2:
        for (int s = 0; s <= nBasis; s++) {
            for (int nx = 0; nx <= s; nx++) {
                for (int ny = 0; ny <= s; ny++) {
                    if(nx + ny == s){
                        for (int spin = 1; spin >= -1; spin -= 2) {
                            state(0) = spin;
                            state(1) = nx;
                            state(2) = ny;
                            states.push_back(state);
                        }
                    }
                }
            }
        }
        break;
    case 3:
        for (int s = 0; s <= nBasis; s++) {
            for (int nx = 0; nx <= s; nx++) {
                for (int ny = 0; ny <= s; ny++) {
                    for (int nz = 0; nz <= s; nz++) {
                        if(nx + ny + nz == s){
                            for (int spin = 1; spin >= -1; spin -= 2) {
                                state(0) = spin;
                                state(1) = nx;
                                state(2) = ny;
                                state(3) = nz;
                                states.push_back(state);
                            }
                        }
                    }
                }
            }
        }
        break;
    }
    cout << states.size() << " orbitals created" << endl;

#ifdef DEBUG
    cout << "void Basis::createCartesianBasis()" << endl;
    for (int i = 0; i < (int)states.size(); i++) {
        for (int j = 0; j < (int)states[0].n_elem; j++) {
            cout << states[i][j] << " ";
        }
        cout << endl;
    }
#endif
}
//------------------------------------------------------------------------------
void Basis::createPolarBasis()
{
    vec state(dim +1);
    switch (dim) {
    case 1:
        cerr << "A polar basis in 1d does not make sense..." << endl;
        exit(1);
        break;
    case 2:
        for (int k = 0; k <= nBasis; k++) {
            for (int i = 0; i <= floor(k / 2); i++) {
                for (int l = -k; l <= k; l++) {
                    if ((2 * i + abs(l) + 1) == k) {

                        // Spin up
                        state(0) = 1;
                        state(1) = i;
                        state(2) = l;
                        states.push_back(state);

                        // Spin down
                        state(0) = -1;
                        state(1) = i;
                        state(2) = l;
                        states.push_back(state);
                    }
                }
            }
        }
        break;
    }
    cout << states.size() << " orbitals created" << endl;
#ifdef DEBUG
    cout << "void Basis::createPolarBasis()" << endl;
    for (int i = 0; i < (int)states.size(); i++) {
        for (int j = 0; j < (int)states[0].n_elem; j++) {
            cout << states[i][j] << " ";
        }
        cout << endl;
    }
#endif
}
//------------------------------------------------------------------------------
const vector<vec> &Basis::getBasis() const
{
    return states;
}
//------------------------------------------------------------------------------
const cx_mat &Basis::getInitalOrbitals() const
{
    return C;
}
//------------------------------------------------------------------------------
const vec &Basis::getX() const
{
    return x;
}
//------------------------------------------------------------------------------
Basis::~Basis()
{
}
//------------------------------------------------------------------------------
