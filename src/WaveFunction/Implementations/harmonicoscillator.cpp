#include "harmonicoscillator.h"

//------------------------------------------------------------------------------
HarmonicOscillator::HarmonicOscillator(Config *cfg, vec quantumNumbers):
    Wavefunction(cfg, quantumNumbers)
{
    try {
        dim = cfg->lookup("system.dim");
        w = cfg->lookup("potential.w");
        aa = cfg->lookup("potential.a");
        aa *= aa;
    } catch (const SettingNotFoundException &nfex) {
        cerr << "OrbitalHarmonicOscillator::Error reading from 'systemSettings' object setting." << endl;
    }
    sqrtW = sqrt(w);
#if DEBUG
//    cout << "OrbitalHarmonicOscillator::OrbitalHarmonicOscillator(Setting* systemSettings, vec quantumNumbers)" << endl;
//    cout << "w = " << w << endl;
//    cout << "dim = " << dim << endl;
#endif
}
//------------------------------------------------------------------------------
double HarmonicOscillator::hermitePolynomial(const int n, const double x)
{
    double hermite_polynomial;

    switch(n){
    case 0:
        hermite_polynomial = 1;
        break;
    case 1:
        hermite_polynomial = 2*x;
        break;
    case 2:
        hermite_polynomial = 4*pow(x,2)  - 2;
        break;
    case 3:
        hermite_polynomial = 8*pow(x,3)  - 12*x;
        break;
    case 4:
        hermite_polynomial = 16*pow(x,4)  - 48*pow(x,2)  + 12;
        break;
    case 5:
        hermite_polynomial = 32*pow(x,5)  - 160*pow(x,3)  + 120*x;
        break;
    case 6:
        hermite_polynomial = 64*pow(x,6)  - 480*pow(x,4)  + 720*pow(x,2)  - 120;
        break;
    case 7:
        hermite_polynomial = 128*pow(x,7)  - 1344*pow(x,5)  + 3360*pow(x,3)  - 1680*x;
        break;
    case 8:
        hermite_polynomial = 256*pow(x,8)  - 3584*pow(x,6)  + 13440*pow(x,4)  - 13440*pow(x,2)  + 1680;
        break;
    case 9:
        hermite_polynomial = 512*pow(x,9)  - 9216*pow(x,7)  + 48384*pow(x,5)  - 80640*pow(x,3)  + 30240*x;
        break;
    case 10:
        hermite_polynomial = 1024*pow(x,10)  - 23040*pow(x,8)  + 161280*pow(x,6)  - 403200*pow(x,4)  + 302400*pow(x,2)  - 30240;
        break;
    case 11:
        hermite_polynomial = 2048*pow(x,11)  - 56320*pow(x,9)  + 506880*pow(x,7)  - 1774080*pow(x,5)  + 2217600*pow(x,3)  - 665280*x;
        break;
    case 12:
        hermite_polynomial = 4096*pow(x,12)  - 135168*pow(x,10)  + 1520640*pow(x,8)  - 7096320*pow(x,6)  + 13305600*pow(x,4)  - 7983360*pow(x,2)  + 665280;
        break;
    case 13:
        hermite_polynomial = 8192*pow(x,13)  - 319488*pow(x,11)  + 4392960*pow(x,9)  - 26357760*pow(x,7)  + 69189120*pow(x,5)  - 69189120*pow(x,3)  + 17297280*x;
        break;
    case 14:
        hermite_polynomial = 16384*pow(x,14)  - 745472*pow(x,12)  + 12300288*pow(x,10)  - 92252160*pow(x,8)  + 322882560*pow(x,6)  - 484323840*pow(x,4)  + 242161920*pow(x,2)  - 17297280;
        break;
    case 15:
        hermite_polynomial = 32768*pow(x,15)  - 1720320*pow(x,13)  + 33546240*pow(x,11)  - 307507200*pow(x,9)  + 1383782400*pow(x,7)  - 2905943040*pow(x,5)  + 2421619200*pow(x,3)  - 518918400*x;
        break;
    case 16:
        hermite_polynomial = 65536*pow(x,16)  - 3932160*pow(x,14)  + 89456640*pow(x,12)  - 984023040*pow(x,10)  + 5535129600*pow(x,8)  - 15498362880*pow(x,6)  + 19372953600*pow(x,4)  - 8302694400*pow(x,2)  + 518918400;
        break;
    case 17:
        hermite_polynomial = 131072*pow(x,17)  - 8912896*pow(x,15)  + 233963520*pow(x,13)  - 3041525760*pow(x,11)  + 20910489600*pow(x,9)  - 75277762560*pow(x,7)  + 131736084480*pow(x,5)  - 94097203200*pow(x,3)  + 17643225600*x;
        break;
    case 18:
        hermite_polynomial = 262144*pow(x,18)  - 20054016*pow(x,16)  + 601620480*pow(x,14)  - 9124577280*pow(x,12)  + 75277762560*pow(x,10)  - 338749931520*pow(x,8)  + 790416506880*pow(x,6)  - 846874828800*pow(x,4)  + 317578060800*pow(x,2)  - 17643225600;
        break;
    case 19:
        hermite_polynomial = 524288*pow(x,19)  - 44826624*pow(x,17)  + 1524105216*pow(x,15)  - 26671841280*pow(x,13)  + 260050452480*pow(x,11)  - 1430277488640*pow(x,9)  + 4290832465920*pow(x,7)  - 6436248698880*pow(x,5)  + 4022655436800*pow(x,3)  - 670442572800*x;
        break;
    case 20:
        hermite_polynomial = 1048576*pow(x,20)  - 99614720*pow(x,18)  + 3810263040*pow(x,16)  - 76205260800*pow(x,14)  + 866834841600*pow(x,12)  - 5721109954560*pow(x,10)  + 21454162329600*pow(x,8)  - 42908324659200*pow(x,6)  + 40226554368000*pow(x,4)  - 13408851456000*pow(x,2)  + 670442572800;
        break;
    default:
        std::cout << "Error. Analytic expressions for the Hermite Polynomial of this degree has not been implemented.";
        exit(1);
        break;
    }

    return hermite_polynomial;
}
//------------------------------------------------------------------------------
