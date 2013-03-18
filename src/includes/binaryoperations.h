#ifndef BINARYOPERATIONS_H
#define BINARYOPERATIONS_H

#include <src/includes/defines.h>
#include <bitset>

using namespace std;

// Binary operations
//void addParticle(const int, bitset<BITS> &state);
//void removeParticle(const int, bitset<BITS> &state);
//int sign(const int, const bitset<BITS> &state);

//------------------------------------------------------------------------------
inline int sign(const int n, const bitset<BITS> &state)
{
    int s = 1;
    for (int i = 0; i < n; i++) {
        if (state[i] != 0) {
            s *= -1;
        }
    }
    return s;
}
//------------------------------------------------------------------------------
inline void addParticle(const int n, bitset<BITS> &state)
{
    bitset<BITS> a;
    a.set(n);
    bitset<BITS> comp = a & state;

    if (comp.count() == false)
        state.set(n);
    else
        state.set(BITS - 1);

}
//------------------------------------------------------------------------------
inline void removeParticle(const int n, bitset<BITS> &state)
{
    bitset<BITS> a;
    a.set(n);
    bitset<BITS> comp = a & state;

    if (comp.count() == true)
        state.set(n, 0);
    else
        state.set(BITS - 1);

}
//------------------------------------------------------------------------------
inline int signBool(const int n, const bool *state)
{
    int s = 1;
    for (int i = 0; i < n; i++) {
        if (state[i]) {
            s *= -1;
        }
    }
    return s;
}
//------------------------------------------------------------------------------
inline void addParticleBool(const int n, bool *state)
{
    if(!state[n])
        state[n] = true;
    else
        state[BITS-1] = true;
}
//------------------------------------------------------------------------------
inline void removeParticleBool(const int n, bool *state)
{
    if(state[n])
        state[n] = false;
    else
        state[BITS-1] = true;
}
//------------------------------------------------------------------------------
#endif // BINARYOPERATIONS_H
