#include "wavefunction.h"
//------------------------------------------------------------------------------
Wavefunction::Wavefunction(Config *cfg, vec quantumNumbers):
    cfg(cfg),
    quantumNumbers(quantumNumbers)
{
    try {
        dim = cfg->lookup("system.dim");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "Wavefunction::Error reading from 'systemSettings' object setting." << endl;
    }
}
//------------------------------------------------------------------------------
Wavefunction::~Wavefunction()
{
}
//------------------------------------------------------------------------------

