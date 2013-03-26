#include "timepropagation.h"

//------------------------------------------------------------------------------
TimePropagation::TimePropagation(Config *cfg):
    cfg(cfg)
{
    try{
        dt = cfg->lookup("timeIntegration.dt");
        N = cfg->lookup("timeIntegration.N");
        cfg->lookupValue("systemSettings.filePath", filePath);
        cfg->lookupValue("systemSettings.filePath", filePath);
        saveToFileInterval = cfg->lookup("systemSettings.saveToFileInterval");
        printProgress  = cfg->lookup("systemSettings.printProgress");
        saveEveryTimeStep = cfg->lookup("systemSettings.saveEveryTimeStep");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "TimePropagation::TimePropagation(Config *cfg)::Error reading parameter from config object." << endl;
        exit(EXIT_FAILURE);
    }
    step = 0;
    t = 0;
    E = zeros(N);
    K = vec(1);
    time = zeros(N);

    filenameOrbitals = filePath + "C.mat";
    filenameSlaterDet = filePath + "A.mat";
    filenameEnergy = filePath + "E.mat";
    filenameDeltaE = filePath + "dE.mat";
    filenameSvdRho = filePath + "svdRho.mat";
    filenameRho = filePath + "rho.mat";
    filenameCorrelation = filePath + "K.mat";
    filenameT = filePath + "time.mat";

#ifdef DEBUG
    cout << "TimePropagation::TimePropagation(Config *cfg)::" << endl
         << "dt \t= " << dt << endl
         << "N \t= " << N << endl;
#endif
}
//------------------------------------------------------------------------------
void TimePropagation::doTimePropagation()
{
    cout << "\n" << endl;
    bool accepted;

    for(step=0; step < N; step++){
        cout << "before step: " << step << " N = " << N << endl;
        accepted = this->stepForward();

        // Saving C and A to disk
        if(accepted){
            time(step) = t;

            // Updating the one-body- and interaction-elements
            V->computeNewElements(C);
            h->computeNewElements(C);

            // Collecting data
            E(step) = slater->getEnergy(A);
            cout << "bla= " << step << endl;
            rho = &orbital->reCalculateRho1(A);
            K = orbital->getCorrelation();
            svdRho = orbital->getSvdRho1();

            saveProgress(step);
            printProgressToScreen(step);
            t += dt;
        }
    }

    // Saving results
    time.save(filenameT, arma_ascii);
}
//------------------------------------------------------------------------------
void TimePropagation::setDependencies(SlaterEquation *slater,
                                            OrbitalEquation *orbital,
                                            Interaction *V,
                                            SingleParticleOperator *h)
{
    this->slater = slater;
    this->orbital = orbital;
    this->V = V;
    this->h = h;
}
//------------------------------------------------------------------------------
void TimePropagation::setInititalState(cx_vec &A, cx_mat &C)
{
    this->A = A;
    this->C = C;
}
//------------------------------------------------------------------------------
cx_mat TimePropagation::getCurrentC()
{
    return C;
}
//------------------------------------------------------------------------------
cx_vec TimePropagation::getCurrentA()
{
    return A;
}
//------------------------------------------------------------------------------
void TimePropagation::printProgressToScreen(uint counter)
{
    if( printProgress ){
        cout << "---------------------------------------------------------------\n";
        cout << "step = " << step << endl
             << "E = " << E(counter) << endl
             << "k = " << K(0) << endl
             << "dt = " << dt << endl
             << "svdRho1 = ";
        svdRho.raw_print();
        cout << endl;
    }
}
//------------------------------------------------------------------------------
void TimePropagation::saveProgress(uint counter)
{
    E.save(filePath + "t_E", arma_ascii);
    if(saveEveryTimeStep){
        stringstream fileName;
        fileName << filePath << "t_C" << counter << ".mat";
        C.save(fileName.str());
        fileName.str("");
        fileName << filePath << "t_A" << counter << ".mat";
        A.save(fileName.str());
        fileName.str("");
        fileName << filePath << "t_rho" << counter << ".mat";
         (*rho).save(fileName.str());
    }
}
//------------------------------------------------------------------------------

