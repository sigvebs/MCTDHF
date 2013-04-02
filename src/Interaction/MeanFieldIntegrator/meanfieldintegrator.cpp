#include "meanfieldintegrator.h"

//------------------------------------------------------------------------------
MeanFieldIntegrator::MeanFieldIntegrator(Config *cfg):
    cfg(cfg)
{
    try{
        nGrid = cfg->lookup("spatialDiscretization.nGrid");
        nOrbitals = cfg->lookup("spatialDiscretization.nSpatialOrbitals");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "MeanFieldIntegrator::MeanFieldIntegrator(Config *cfg)"
             << "::Error reading from config object." << endl;
    }

    V2 = field<cx_vec>(nOrbitals, nOrbitals);

    for(int p=0; p < nOrbitals; p++){
        for(int q=0; q < nOrbitals; q++){
            V2(p,q) = zeros<cx_vec>(nGrid);
        }
    }

    // MPI-------------------------------------------
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &nNodes);

    int tot = 0.5*nOrbitals*(nOrbitals + 1);

    allQR = imat(nOrbitals, nOrbitals);
    int node = 0;
    int s = 0;
    for (int q = 0; q < nOrbitals; q++) {
        for (int r = q; r < nOrbitals; r++) {
            if (myRank == node){
                myQR.push_back(pair<int, int>(q,r));
            }
            allQR(q,r) = node;
            s++;
            if(s >= BLOCK_SIZE(node, nNodes, tot)){
                s = 0;
                node++;
            }
        }
    }


    MPI_Barrier(MPI_COMM_WORLD);
    for(pair<int,int>qr: myQR){
        cout << myRank << " " << qr.first << " " << qr.second << endl;
    }

    if(myRank == 0){
        cout << "nOrbitals = " << nOrbitals << endl;
        cout << "tot = " << tot << endl;
        cout << allQR << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
//    exit(1);
}
//------------------------------------------------------------------------------
void MeanFieldIntegrator::computeMeanField(const cx_mat &C)
{
    int q, r;
    for(pair<int,int>qr: myQR){
        q = qr.first;
        r = qr.second;
        integrate(q, r, C, V2(q, r));
    }

    for (int q = 0; q < nOrbitals; q++) {
        for (int r = q; r < nOrbitals; r++) {
            MPI_Bcast( V2(q,r).memptr(), nGrid , MPI_DOUBLE_COMPLEX, allQR(q,r), MPI_COMM_WORLD ) ;
            V2(r, q) = conj(V2(q, r));
        }
    }

//    if(myRank == 1){
//        for (int q = 0; q < nOrbitals; q++) {
//            for (int r = q; r < nOrbitals; r++) {
//                cout << V2(q,r) << endl;
//            }
//        }
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    exit(1);

//    MPI_Barrier(MPI_COMM_WORLD);
//    exit(1);

//    for (int q = 0; q < nOrbitals; q++) {
//        integrate(q, q, C, V2(q,q));

//        for (int r = q+1; r < nOrbitals; r++) {
//            integrate(q, r, C, V2(q,r));
//            V2(r,q) = conj(V2(q,r));
//        }
    //    }
}

//------------------------------------------------------------------------------
void MeanFieldIntegrator::cleanUp()
{
    while (!potential.empty())
    {
      InteractionPotential* pot = potential.back();
      potential.pop_back();
      delete pot;
    }
}
//------------------------------------------------------------------------------
void MeanFieldIntegrator::addPotential(InteractionPotential *interactionPot)
{
    potential.push_back(interactionPot);
}
//------------------------------------------------------------------------------
MeanFieldIntegrator::~MeanFieldIntegrator()
{

}
//------------------------------------------------------------------------------
