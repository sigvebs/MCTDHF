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
        saveToFileInterval = cfg->lookup("systemSettings.saveToFileIntervalTime");
        printProgress  = cfg->lookup("systemSettings.printProgress");
        saveEveryTimeStep = cfg->lookup("systemSettings.saveEveryTimeStep");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "TimePropagation::TimePropagation(Config *cfg)::Error reading parameter from config object." << endl;
        exit(EXIT_FAILURE);
    }

    step = 0;
    t = 0;
    E = zeros(N/saveToFileInterval+1);
    K = vec(1);
    time = zeros(N);

    filenameOrbitals = filePath + "/C.mat";
    filenameSlaterDet = filePath + "/A.mat";
    filenameEnergy = filePath + "/E.mat";
    filenameDeltaE = filePath + "/dE.mat";
    filenameSvdRho = filePath + "/svdRho.mat";
    filenameRho = filePath + "/rho.mat";
    filenameCorrelation = filePath + "/K.mat";
    filenameT = filePath + "/time.mat";

    isMaster = true;
#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &nNodes);

    isMaster = (bool)(myRank == 0);
#endif
#ifdef DEBUG
    cout << "TimePropagation::TimePropagation(Config *cfg)::" << endl
         << "dt \t= " << dt << endl
         << "N \t= " << N << endl;
#endif
}
//------------------------------------------------------------------------------
void TimePropagation::doTimePropagation()
{
    int counter = 0;
    bool accepted;

    for(step=0; step < N; step++){
        accepted = this->stepForward();

        // Saving C and A to disk
        if(accepted){
            if((step % saveToFileInterval == 0 || step == N-1) ){
                time(counter) = t;

#ifdef USE_MPI
                MPI_Bcast( C.memptr(), C.n_elem , MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD );
#endif

                // Updating the one-body- and interaction-elements
                V->computeNewElements(C);
                h->computeNewElements(C, t);

                // Collecting data
                if(isMaster){
                    E(counter) = slater->getEnergy(A);
                    rho = &orbital->reCalculateRho1(A);
                    K = orbital->getCorrelation();
                    svdRho = orbital->getSvdRho1();

                    saveProgress(counter);
                    printProgressToScreen(counter);
                }
                counter++;
            }
            t += dt;
        }
    }

    // Saving results
    if(isMaster)
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
void TimePropagation::renormalize(cx_mat &D)
{
    // Re-normalization of C using SVD
    cx_mat X;
    vec s;
    cx_mat Y;
    svd_econ(X, s, Y, D);
    D = X*Y.t();
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
             << "t = " << t << endl
             << "dt = " << dt << endl
             << "svdRho1 = ";
        svdRho.raw_print();
        cout << endl;
    }
}
//------------------------------------------------------------------------------
void TimePropagation::saveProgress(uint counter)
{
    E.save(filePath + "/t_E.mat", arma_ascii);
    if(saveEveryTimeStep){
        stringstream fileName;
        fileName << filePath << "/t_C" << counter << ".mat";
        C.save(fileName.str());
        fileName.str("");
        fileName << filePath << "/t_A" << counter << ".mat";
        A.save(fileName.str());
        fileName.str("");
        fileName << filePath << "/t_rho" << counter << ".mat";
         (*rho).save(fileName.str());
    }
}
//------------------------------------------------------------------------------
