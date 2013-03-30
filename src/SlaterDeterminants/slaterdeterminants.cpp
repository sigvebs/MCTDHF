#include "slaterdeterminants.h"

//------------------------------------------------------------------------------
SlaterDeterminants::SlaterDeterminants(Config *cfg, vector<vec> sps): cfg(cfg), sps(sps)
{
    try {
        nParticles = cfg->lookup("system.nParticles");
        conservedEigenSpin = cfg->lookup("system.conserveSpin");
        conservedEigenSpinValue = cfg->lookup("system.spinValue");
        cfg->lookupValue("systemSettings.filePath", filePath);
    } catch (const SettingNotFoundException &nfex) {
        cerr << "SlaterDeterminants(Config *cfg, vector<vec> sps)::Error reading config object." << endl;
        exit(EXIT_FAILURE);
    }
#ifdef DEBUG
    cout << "SlaterDeterminants(Config *cfg, vector<vec> sps)" << endl;
    cout << "nParticles = " << nParticles << endl;
    cout << "conservedEigenSpin = " << conservedEigenSpin << endl;
    cout << "conservedEigenSpinValue = " << conservedEigenSpinValue << endl;
#endif
}
//------------------------------------------------------------------------------
void SlaterDeterminants::createSlaterDeterminants()
{
    // Creating all possible states
    int nSps = sps.size();
    vec state = zeros(nParticles);

    // Creating an initial state
    for (int i = 0; i < nParticles; i++) {
        state(i) = i;
    }

    if (!conservedEigenSpin || checkEigenSpin(state)) {
        binStates.push_back(createBinaryState(state));
    }

    // Creating all possible slater determinants.
    while (true) {
        state = odometer(state, nSps, nParticles);

        if (sum(state) == 0)
            break;

        // Storing accepted states.
        if (!conservedEigenSpin || checkEigenSpin(state)) {
            binStates.push_back(createBinaryState(state));
        }
    }

    cout << binStates.size()  << " Slater determinants created." << endl;

    saveSlaterDeterminantsToDisk();
#ifdef DEBUG
    cout << "SlaterDeterminants::createSlaterDeterminants()" << endl;
    for(int i=0; i<(int)binStates.size(); i++){
        cout << binStates[i] << endl;
    }
#endif
}
//------------------------------------------------------------------------------
bool SlaterDeterminants::checkEigenSpin(vec state)
{
    bool eigenSpin = true;
    int sumEigenSpin = 0;

    // Summing up the total eigenspin
    for (int i = 0; i < nParticles; i++) {
        sumEigenSpin += (sps[state(i)])[0];
    }

    // Checking if the total eigenspin of a state is conserved
    if (sumEigenSpin != conservedEigenSpinValue) {
        eigenSpin = false;
    }

    return eigenSpin;
}
//------------------------------------------------------------------------------
bitset<BITS> SlaterDeterminants::createBinaryState(vec state)
{
    // Constructing the binary representation
    bitset<BITS> binState;

    for (int i = 0; i < (int)state.size(); i++) {
        try {
            binState.set(state(i));
        } catch (exception& e) {
            cout << "Exception: " << e.what() << endl;
            cout << "To run with this basis the number of bits must be increased. Current number of bits are " << BITS
                 << " , increase BITS > " << state(i) - 1 << ". Change BITS in defines.h and recompile." << endl;
            exit(1);
        }
    }

    return binState;
}
//------------------------------------------------------------------------------
vec SlaterDeterminants::odometer(const vec &oldState, int Nsp, int N)
{
    vec newState = oldState;
    double l;

    for (int j = N - 1; j >= 0; j--) {
        if (newState(j) < Nsp - N + j) {
            l = newState(j);
            for (int k = j; k < N; k++) {
                newState(k) = l + 1 + k - j;
            }
            return newState;
        }
    }

    newState = zeros(nParticles);
    return newState;
}
//------------------------------------------------------------------------------
const vector<bitset<BITS> > &SlaterDeterminants::getSlaterDeterminants() const
{
    return binStates;
}
//------------------------------------------------------------------------------
void SlaterDeterminants::saveSlaterDeterminantsToDisk()
{
    ofstream SdFile;
    SdFile.open (filePath + "/SlaterDeterminants.sd");

    SdFile << (int)binStates.size() << " " << BITS << endl;
    for(int i=0; i<(int)binStates.size(); i++){
        SdFile << binStates[i] << endl;
    }
    SdFile.close();
}
//------------------------------------------------------------------------------
cx_vec SlaterDeterminants::getCoefficients()
{
    return A;
}
//------------------------------------------------------------------------------
void SlaterDeterminants::createInitialState()
{
    A = randu<cx_vec>(binStates.size());
    A = A/sqrt(cdot(A, A));
}
//------------------------------------------------------------------------------
void SlaterDeterminants::load()
{
    string loadPath;
    string filenameA;
    string filenameSlatderDet;
    try{
        cfg->lookupValue("loadDataset.loadDatasetPath", loadPath);
        cfg->lookupValue("loadDataset.A", filenameA);
        cfg->lookupValue("loadDataset.slatderDet", filenameSlatderDet);
    }catch (const SettingNotFoundException &fioex) {
        cerr << "Basis::loadOrbitals()::Loadpath not found in config file." << endl;
        exit(EXIT_FAILURE);
    }

    // Laoding the Slater determinants
    int nSlater;
    int bits;
    vec state = zeros(nParticles);
    ifstream SdFile(loadPath + "/" + filenameSlatderDet);
    string line;
    if (SdFile.is_open())
    {
        // Assuming the Slater determinant file has the structure
        // nSlater BITS
        // 001110000000000000000000000
        // 001100100000000000000000000
        // 000001001000000000000000000
        // 00...

        getline (SdFile, line);
        string word;
        stringstream stream(line);
        getline(stream, word, ' ');
        nSlater = atoi(word.c_str());
        getline(stream, word, ' ');
        bits = atoi(word.c_str());

        for(int n=0; n<nSlater; n++){
            state.zeros();
            getline (SdFile, line);

            int counter = 0;
            for(int i= bits-1; i>=0; i--){
                int occupied = line[i] - 48;

                if(occupied){
                    state(counter++) = bits - i - 1;
                }
            }
            binStates.push_back(createBinaryState(state));
        }
        SdFile.close();
    }else{
        cerr << "SlaterDeterminants::load()::Error reading Slater determint from file." << endl;
        exit(EXIT_FAILURE);
    }


    // Loading the coefficients
    A.load(loadPath + "/" + filenameA);
}
//------------------------------------------------------------------------------
