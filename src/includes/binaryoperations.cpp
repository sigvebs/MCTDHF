#include "binaryoperations.h"
//------------------------------------------------------------------------------
int sign(const int n, const bitset<BITS> &state)
{
    if (!state[BITS - 1]) {
        int s = 1;
        for (int i = 0; i < n; i++) {
            if (state[i] != 0) {
                s *= -1;
            }
        }
        return s;
    }
    return 1;
}
//------------------------------------------------------------------------------
void addParticle(const int n, bitset<BITS> &state)
{
    if (!state[BITS - 1]) {
        bitset<BITS> a;
        a.set(n);
        bitset<BITS> comp = a & state;

        if (comp.count() == false)
            state.set(n);
        else
            state.set(BITS - 1);
    }
}
//------------------------------------------------------------------------------
void removeParticle(const int n, bitset<BITS> &state)
{
    if (!state[BITS - 1]) {
        bitset<BITS> a;
        a.set(n);
        bitset<BITS> comp = a & state;

        if (comp.count() == true)
            state.set(n, 0);
        else
            state.set(BITS - 1);
    }
}
//------------------------------------------------------------------------------
