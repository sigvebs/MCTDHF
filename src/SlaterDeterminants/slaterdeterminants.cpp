#include "slaterdeterminants.h"

//------------------------------------------------------------------------------
SlaterDeterminants::SlaterDeterminants(Config *cfg, vector<vec> sps): cfg(cfg), sps(sps)
{ try {
        nParticles = cfg->lookup("system.nParticles");
        conservedEigenSpin = cfg->lookup("system.conserveSpin");
        conservedEigenSpinValue = cfg->lookup("system.spinValue");
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
        binStatesBool.push_back(createBinaryStateBool(state));
        binStates.push_back(createBinaryState(state));
    }

    // Creating all possible slater determinants.
    while (true) {
        state = odometer(state, nSps, nParticles);

        if (sum(state) == 0)
            break;

        // Storing accepted states.
        if (!conservedEigenSpin || checkEigenSpin(state)) {
            binStatesBool.push_back(createBinaryStateBool(state));
            binStates.push_back(createBinaryState(state));
        }
    }

    cout << binStates.size()  << " Slater determinants created." << endl;

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
bool *SlaterDeterminants::createBinaryStateBool(vec state)
{
    // Constructing the binary representation
    bool *binState = new bool[BITS];
    for(int i=0; i<BITS; i++)
        binState[i] = 0;

    for (int i = 0; i < (int)state.size(); i++) {
        try {
            binState[(int)state(i)] = 1;
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
const vector<bool *> &SlaterDeterminants::getSlaterDeterminantsBool() const
{
    return binStatesBool;
}
//------------------------------------------------------------------------------
