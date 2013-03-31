#include "grid.h"

//------------------------------------------------------------------------------
Grid::Grid(Config *cfg):
    cfg(cfg)
{
    try{
        dim = cfg->lookup("system.dim");
        cfg->lookupValue("systemSettings.filePath", path);
        periodicBoundaries=cfg->lookup("spatialDiscretization.periodicBoundaries");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "Grid(Config *cfg)::Error reading from config object." << endl;
        exit(EXIT_FAILURE);
    }
    nGrid = 0;
    nGridX = 0;
    nGridY = 0;
    nGridZ = 0;
    DX = 0;
    DY = 0;
    DZ = 0;
}
//------------------------------------------------------------------------------
void Grid::createInitalDiscretization()
{
    switch (dim) {
    case 1:
        oneDimDiscretization();
        break;
    case 2:
        twoDimDiscretization();
        break;
    case 3:
        threeDimDiscretization();
        break;
    }

    Setting &root = cfg->getRoot();
    Setting &tmp = root["spatialDiscretization"];
    tmp.remove("nGrid");
    tmp.add("nGrid", Setting::TypeInt) = (int)R.n_cols;
}
//------------------------------------------------------------------------------
void Grid::oneDimDiscretization()
{
    vec X;
    double Lx;
    try{
        Lx = cfg->lookup("spatialDiscretization.oneD.Lx");
        nGridX = cfg->lookup("spatialDiscretization.oneD.nGridX");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "Grid::oneDimDiscretization()::Error reading from config object." << endl;
        exit(EXIT_FAILURE);
    }

    nGrid = nGridX;
    R = zeros(dim, nGrid);

    if(periodicBoundaries){
        DX = 2.0*Lx/(double)(nGridX);
        X = linspace<vec>(-Lx, Lx - DX,nGridX);
    }else{
        X = linspace<vec>(-Lx, Lx, nGridX);
        DX = X(1) - X(0);
    }

    for(int i=0; i<nGrid; i++){
        R(0, i) = X(i);
    }
}
//------------------------------------------------------------------------------
void Grid::twoDimDiscretization()
{
    double Lx, Ly;
    try{
        Lx = cfg->lookup("spatialDiscretization.twoD.Lx");
        Ly = cfg->lookup("spatialDiscretization.twoD.Ly");
        nGridX = cfg->lookup("spatialDiscretization.twoD.nGridX");
        nGridY = cfg->lookup("spatialDiscretization.twoD.nGridY");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "Grid::twoDimDiscretization()::Error reading from config object." << endl;
        exit(EXIT_FAILURE);
    }

    nGrid = nGridX*nGridY;

    R = zeros(dim, nGrid);

    vec X, Y;

    if(periodicBoundaries){
        DX = 2.0*Lx/(double)(nGridX);
        DY = 2.0*Ly/(double)(nGridY);
        X = linspace<vec>(-Lx, Lx - DX, nGridX);
        Y = linspace<vec>(-Lx, Lx - DY, nGridY);
    }else{
        X = linspace<vec>(-Ly, Ly, nGridX);
        Y = linspace<vec>(-Ly, Ly, nGridY);
        DX = X(1) - X(0);
        DY = Y(1) - Y(0);
    }

    int counter = 0;
    for(int i=0; i<nGridX; i++){
        for(int j=0; j<nGridY; j++){
            vec2 r;
            r(0) = X(i);
            r(1) = Y(j);
            R.col(counter++) = r;
        }
    }
}
//------------------------------------------------------------------------------
void Grid::threeDimDiscretization()
{
    cerr << "Grid::threeDimDiscretization():";
    cerr << "3D not yet implemented!";
    exit(EXIT_FAILURE);
}
//------------------------------------------------------------------------------
const vec Grid::at(int i) const
{
    return R.col(i);
}
//------------------------------------------------------------------------------
void Grid::saveGrid()
{
    R.save(path + "/r.mat");
}
//------------------------------------------------------------------------------
void Grid::loadGrid()
{
    string loadPath;
    string filenameR;
    try{
        cfg->lookupValue("loadDataset.loadDatasetPath", loadPath);
        cfg->lookupValue("loadDataset.R", filenameR);
    }catch (const SettingNotFoundException &fioex) {
        cerr << "Basis::loadOrbitals(string path)::Loadpath not found in config file." << endl;
        exit(EXIT_FAILURE);
    }

    R.load(filenameR);

    switch (dim) {
    case 1:
        DX = R(0,1) - R(0,0);
        break;
    case 2:
        DX = R(0,1) - R(0,0);
        DY = R(1,1) - R(1,0);
        break;
    case 3:
        DX = R(0,1) - R(0,0);
        DY = R(1,1) - R(1,0);
        DZ = R(2,1) - R(2,0);
        break;
    }

    Setting &root = cfg->getRoot();
    Setting &tmp = root["spatialDiscretization"];
    tmp.remove("nGrid");
    tmp.add("nGrid", Setting::TypeInt) = (int)R.n_elem;
}
//------------------------------------------------------------------------------
