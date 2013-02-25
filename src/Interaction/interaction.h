#ifndef INTERACTION_H
#define INTERACTION_H

// Local includes
#include <src/includes/defines.h>
#include <src/includes/lib.h>

// Library includes
#include <armadillo>
#include <vector>
#include <unordered_map>
#include <libconfig.h++>

using namespace libconfig;
using namespace std;
using namespace arma;

class Interaction
{
public:
    Interaction(Config* cfg, vector<vec> orbitals);
    void computeNewElements(const cx_mat &C);
    const cx_double at(const int p, const int q, const int r, const int s);
    const cx_vec &meanField(const int p, const int q);
protected:
    void computeMeanField();
    void computeInteractionelements();
    cx_vec integrate(const int q, const int r);
    cx_double integrate(const int p, const int q, const int r, const int s);

    Config* cfg;

    vector<vec> orbitals;
    vec x;
    double dx;

    int nGrid;
    int nOrbitals;
    int nSpatialOrbitals;

    double aa;
    cx_mat C;

    // Local storage of data
    unordered_map<int, cx_double> interactionElements;
    field<cx_vec> V2;
};

#endif // INTERACTION_H
