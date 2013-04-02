#ifndef DEFINES_H
#define DEFINES_H

#define TESTING 0
#define BITS 31
#define BUF_SIZE 33
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164

enum complexTimeIntegrator_t {
    CT_RUNGE_KUTTA4, CT_RUNGE_KUTTA_FEHLBERG
};
enum timeIntegrator_t {
    RUNGE_KUTTA4, RUNGE_KUTTA_FEHLBERG
};
enum sIntegrator_t {
    MONTE_CARLO, GAUSS_LAGUERRE, GAUSS_HERMITE, INTERACTION_INTEGRATOR, MONTE_CARLO_IS
};
enum coordinateTypes_t {
    CARTESIAN, POLAR
};
enum differentialOpertor {
    DO_FINITE_DIFFERENCE_1d, DO_FINITE_DIFFERENCE_FIVE_POINT_1D, DO_FOURIER_1D,
    DO_FINITE_DIFFERENCE_2d, DO_FiniteDifferenceFivePoint2d
};
enum orbitalBasisType_t{
    OBT_HARMONIC_OSCILLATOR, OBT_HYDROGEN_LIKE, OBT_RAND_UNITARY_MATRIX
};
enum interactionPotential_t{
    IP_HARMONIC_OSCILLATOR, IP_SHIELDED_COULOMB
};
enum meanFieldIntegrator_t{
    MF_TRAPEZODIAL, MF_LOW_RANK_APPROXIMATION
};
enum oneBodyPotential_t{
    HARMONIC_OSCILLATOR_ONE_BODY, COULOMB_INTERACTION_NUCLEUS, SIMPLE_LASER, ANHARMONIC_DOUBLE_WELL
};

#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(i,p,n) (((p)*((i)+1)-1)/(n))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define CEILING(i,j) (((i)+(j)-1)/(j))

#endif // DEFINES_H
