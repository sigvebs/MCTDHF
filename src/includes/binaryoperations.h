#ifndef BINARYOPERATIONS_H
#define BINARYOPERATIONS_H

#include <src/includes/defines.h>
#include <bitset>

using namespace std;

// Binary operations
void addParticle(const int, bitset<BITS> &state);
void removeParticle(const int, bitset<BITS> &state);
int sign(const int, const bitset<BITS> &state);

#endif // BINARYOPERATIONS_H
