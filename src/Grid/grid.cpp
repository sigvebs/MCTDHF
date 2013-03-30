#include "grid.h"

//------------------------------------------------------------------------------
Grid::Grid(Config *cfg):
    cfg(cfg)
{
    try{
        dim = cfg->lookup("system.dim");
        cfg->lookupValue("systemSettings.filePath", path);
        nGrid = cfg->lookup("spatialDiscretization.nGrid");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "Grid(Config *cfg)::Error reading from config object." << endl;
        exit(EXIT_FAILURE);
    }
}
//------------------------------------------------------------------------------
void Grid::createInitalDiscretization()
{
    double L;
    double Lx, Ly, Lz;
    bool periodicBoundaries;
    try{
        L = cfg->lookup("spatialDiscretization.latticeRange");
        periodicBoundaries=cfg->lookup("spatialDiscretization.periodicBoundaries");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "discretizeBasis()::Error reading from config object." << endl;
        exit(EXIT_FAILURE);
    }
    Lx = Ly = Lz = L;
    // Adding dx to the config file
    Setting &root = cfg->getRoot();
    Setting &tmp = root["spatialDiscretization"];

    switch (dim) {
    case 1:
        if(periodicBoundaries){
            DX = 2.0*L/(double)(nGrid);
            X = linspace<vec>(-Lx, Lx - DX,nGrid);
        }else{
            X = linspace<vec>(-Lx, Lx, nGrid);
            DX = X(1) - X(0);
        }
        break;
    case 2:
        if(periodicBoundaries){
            DX = 2.0*L/(double)(nGrid);
            DY = 2.0*L/(double)(nGrid);
            X = linspace<vec>(-Lx, Lx - DX, nGrid);
            Y = linspace<vec>(-Lx, Lx - DY, nGrid);
        }else{
            X = linspace<vec>(-Ly, Ly, nGrid);
            Y = linspace<vec>(-Ly, Ly, nGrid);
            DX = X(1) - X(0);
            DY = Y(1) - Y(0);
        }
        break;
    case 3:
        if(periodicBoundaries){
            DX = 2.0*L/(double)(nGrid);
            DY = 2.0*L/(double)(nGrid);
            DZ = 2.0*L/(double)(nGrid);
            X = linspace<vec>(-Lx, Lx - DX, nGrid);
            Y = linspace<vec>(-Ly, Ly - DY, nGrid);
            Z = linspace<vec>(-Lz, Lz - DZ, nGrid);
        }else{
            X = linspace<vec>(-Lx, Lx, nGrid);
            Y = linspace<vec>(-Ly, Ly, nGrid);
            Z = linspace<vec>(-Lz, Lz, nGrid);
            DX = X(1) - X(0);
            DY = Y(1) - Y(0);
            DZ = Z(1) - Z(0);
        }
        break;
    }

    tmp.add("dx", Setting::TypeFloat) = DX;
    tmp.add("dy", Setting::TypeFloat) = DY;
    tmp.add("dz", Setting::TypeFloat) = DZ;
}
//------------------------------------------------------------------------------
void Grid::saveGrid()
{
    switch (dim) {
    case 1:
        X.save(path + "/x.mat");
        break;
    case 2:
        X.save(path + "/x.mat");
        Y.save(path + "/y.mat");
        break;
    case 3:
        X.save(path + "/x.mat");
        Y.save(path + "/y.mat");
        Z.save(path + "/z.mat");
        break;
    }
}
//------------------------------------------------------------------------------
void Grid::loadGrid()
{
    string loadPath;
    string filenameX, filenameY, filenameZ;
    try{
        cfg->lookupValue("loadDataset.loadDatasetPath", loadPath);
        cfg->lookupValue("loadDataset.X", filenameX);
        cfg->lookupValue("loadDataset.Y", filenameY);
        cfg->lookupValue("loadDataset.Z", filenameZ);
    }catch (const SettingNotFoundException &fioex) {
        cerr << "Basis::loadOrbitals(string path)::Loadpath not found in config file." << endl;
        exit(EXIT_FAILURE);
    }

    switch (dim) {
    case 1:
        X.load(loadPath + "/" + filenameX);
        DX = X(1) - X(0);
        break;
    case 2:
        X.load(loadPath + "/" + filenameX);
        Y.load(loadPath + "/" + filenameY);
        DX = X(1) - X(0);
        DY = Y(1) - Y(0);
        break;
    case 3:
        X.load(loadPath + "/" + filenameX);
        Y.load(loadPath + "/" + filenameY);
        Z.load(loadPath + "/" + filenameZ);
        DX = X(1) - X(0);
        DY = Y(1) - Y(0);
        DZ = Z(1) - Z(0);
        break;
    }

    Setting &root = cfg->getRoot();
    Setting &tmp = root["spatialDiscretization"];
    tmp.add("dx", Setting::TypeFloat) = DX;
    tmp.remove("nGrid");
    tmp.add("nGrid", Setting::TypeInt) = (int)X.n_rows;
}
//------------------------------------------------------------------------------
