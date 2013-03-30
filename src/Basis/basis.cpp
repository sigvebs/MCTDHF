#include "basis.h"
//------------------------------------------------------------------------------
Basis::Basis(Config *cfg):
    cfg(cfg)
{
    string filePath;
    try{
        dim = cfg->lookup("system.dim");
        coordinateType = cfg->lookup("system.coordinateType");
        cfg->lookupValue("systemSettings.filePath", filePath);
        nGrid = cfg->lookup("spatialDiscretization.nGrid");
        nBasis = cfg->lookup("system.shells");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "Basis(Config *cfg)::Error reading from config object." << endl;
        exit(EXIT_FAILURE);
    }

    nSpatialOrbitals = states.size()/2;

#ifdef DEBUG
    cout << "Basis::Basis(Config *cfg)" << endl
         << "coordinateType \t= " << coordinateType << endl
         << "dim \t\t= " << dim << endl;
#endif
}
//------------------------------------------------------------------------------
void Basis::createBasis()
{
    states.clear();

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
    try{
        tmp.remove("nSpatialOrbitals");
    }catch(const SettingNotFoundException &nfex) {
    }
    tmp.add("nSpatialOrbitals", Setting::TypeInt) = nSpatialOrbitals;

    cout << "--------------------_" << endl;
    for(vec state:states){
        cout << state << endl;
    }
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
void Basis::SVD(cx_mat &D)
{
    // Forcing orthogonality by performing a SVD decomposition
    cx_mat X;
    vec s;
    cx_mat Y;
    svd_econ(X, s, Y, D);
    D = X*Y.t();
}
//------------------------------------------------------------------------------
const vector<vec> &Basis::getBasis() const
{
    return states;
}
//------------------------------------------------------------------------------
const field<cx_mat> &Basis::getOrbitals() const
{
    return C;
}
//------------------------------------------------------------------------------
void Basis::setGrid(Grid *grid)
{
    this->grid = grid;
}
//------------------------------------------------------------------------------
void Basis::loadOrbitals()
{

    string loadPath;
    string filenameC;
    try{
        cfg->lookupValue("loadDataset.loadDatasetPath", loadPath);
        cfg->lookupValue("loadDataset.C", filenameC);
    }catch (const SettingNotFoundException &fioex) {
        cerr << "Basis::loadOrbitals(string path)::Loadpath not found in config file." << endl;
        exit(EXIT_FAILURE);
    }

    Setting &root = cfg->getRoot();
    Setting &tmp2 = root["system"];

    switch (dim) {
    case 1:
        C = field<cx_mat>(1);
        C(0).load(loadPath + "/" + filenameC);

        tmp2.remove("shells");
        tmp2.add("shells", Setting::TypeInt) = (int)C(0).n_cols - 1;
        nBasis =  (int)C(0).n_cols - 1;
        break;
    case 2:
        C = field<cx_mat>(nSpatialOrbitals);
        for(int i=0; i<nSpatialOrbitals; i++){
            stringstream fileName;
            fileName << loadPath << "C" << i << ".mat";
            C(i).load(loadPath + "/" + filenameC);
        }
    }

    // Recreating the orbital states.
    this->createBasis();
}
//------------------------------------------------------------------------------
Basis::~Basis()
{
}
//------------------------------------------------------------------------------
