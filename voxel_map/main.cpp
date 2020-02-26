#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <cmath>
#include <string>
//#include <cstdio>
#include<stdio.h>
#include<sstream>
//#include <ostream> // Don't need this?
// This is buggy. Should I not use qt?

using namespace std;
using std::ofstream; using std::string;

vector<int> int_linspace(int framefirst, int framelast, int Nframes);
vector<vector<vector<double>>> voxellate_basic(int frame, int Nx, int Ny, int Nz, vector<double> x_centres, vector<double> y_centres, vector<double> z_centres, vector<vector<double>> xes, vector<vector<double>> ys, vector<vector<double>> zs); // voxellate(Nx, Ny, Nz, x_centres, y_centres, z_centres, posvecs)

//vector<vector<double>> voxellate_basic(int Nx, int Ny, int Nz, vector<double> x_centres, vector<double> y_centres, vector<double> z_centres, vector<vector<double>> xes, vector<vector<double>> ys, vector<vector<double>> zs); // voxellate(Nx, Ny, Nz, x_centres, y_centres, z_centres, posvecs)

int main()
{
    int M, N, gridspacing;
    double wall, charge, T, Kangle, Kbond;
    string filename;

    int Nvoxels    = 55;
    int totframes  = 10000000/10000;

    // File selection
    //string
    //filename = "chaingrid_quadratic_M9N101_gridspacing40_Langevin_wall1.042_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debyecutoff3_charge-1_T3_theta0is180_twofirst_are_fixed.lammpstrj";
    //filename = "chaingrid_quadratic_M9N101_gridspacing40_Langevin_wall1.042_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_T3_theta0is180_twofirst_are_fixed.lammpstrj";
    //filename = "C:/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/chaingrid_quadratic_M9N101_gridspacing40_Langevin_wall1.042_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_T3_theta0is180_twofirst_are_fixed.lammpstrj";

    // File name for testing (look more into which files I've deleted and which I haven't later):
    filename = "/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/chaingrid_quadratic_M9N101_ljdebye1.042_angle_Langevin_wall1.042_Kangle25_Kbond2000_T310_theta0is180_cutoffs3_sigma1_firstatomfixed.lammpstrj";

    M = 9;
    N = 101;
    wall = 1.042;
    Kangle = 14.0186574854529;
    Kbond = 140.186574854529;
    charge = -1;
    T = 3;
    gridspacing = 40;

    // File output selection:
    int Nout = 3;
    vector<int> outframes = vector<int>(Nout);
    outframes[0] = 0;
    outframes[1] = totframes/2;
    outframes[2] = totframes-1;

    // Readying voxel centres:
    vector<double> ns = vector<double>(Nvoxels);
    for(int i=0; i<Nvoxels; i++){ns[i]=i;}

    ///*
    ifstream inFile;
    inFile.open(filename);
    if (!inFile)
    {
        cerr << "Could not open file :(" << endl;
        exit(1);
    }
    //*/
    //int totframes  = 10000000/10000;
    int framefirst = 0;
    int framelast  = totframes;
    int Nall       = M*N;
    int skiplines  = 9;             // If we hit 'ITEM:', skip this many steps...
    vector<int> indices = vector<int>(totframes);
    //indices = int_linspace(framefirst,framelast,totframes);
    //*
    int framespacing = (framelast-framefirst)/(totframes-1); //
    for (int i=0; i<totframes-1; i++)
    {
        indices[i] = framefirst + i*framespacing;
        cout << "i=" << i << "; value =" << framefirst + i*framespacing << endl;
    }

    // ------------------- Read in the data and store in matrices ------------------- //
    string currline;
    int id, type, mol;
    double x,y,z, vx, vy,vz;
    int miss       = 0;
    int Ncount     = 0;
    int storeframe = 0;
    int framecount = 0;  // To put values into matrix
    int timestep   = -1; // A counter for the time step.
    double minx     = 100;
    double maxx     = -100;
    double miny     = 100;
    double maxy     = -100;
    double minz     = 0;   // We know this
    double maxz     = 0;   // We will find this
    vector<double> maxzs = vector<double>(M);
    /* // Both these formulations work
    double xes[Nout][Nall];
    double ys[Nout][Nall];
    double zs[Nout][Nall];
    */
    vector<vector<double>> xes(Nout, vector<double>(Nall));
    vector<vector<double>> ys(Nout, vector<double>(Nall));
    vector<vector<double>> zs(Nout, vector<double>(Nall));

    while (std::getline(inFile, currline))
    {
        // Extract data
        std::istringstream iss(currline);

        if (!(iss >> id >> type >> mol >> x >> y >> z >> vx >> vy >> vz)) {
            miss = 1;
            cout << "timestep:" << timestep << endl;
            continue;} // if not the line we want
        if (miss==1){
            miss = 0;
            Ncount = 0;
            timestep++; // Update the time step
            storeframe = outframes[framecount]; // Test if this is the frame we want
        }


        // For finding the bounds of the voxel matrix // Should I have this here? Takes more time, but is more robust this way.
        if (z>maxzs[mol-1]){maxzs[mol-1]=z;} // Getting the max value of z for each chain (regardless of the time step. We want the same cutting for all times) -- If we put this below the if-test, this complicates running this file several times for different output times
        if (x>maxx){maxx=x;}
        if (x<minx){minx=x;}
        if (y>maxy){maxy=y;}
        if (y<miny){miny=y;}

        // Process data

        if (!(storeframe==timestep)){continue;} // Don't do anything else if we don't want to store this frame

        xes[framecount][Ncount] = x;
        ys[framecount][Ncount]  = y;
        zs[framecount][Ncount]  = z;

        Ncount++;
        if (Ncount==Nall){
            framecount++;
        }
    }

    // ------------------- Perform the cutting ------------------- //
    // Finding value maxz
    for(int i=0; i<M; i++){maxz+=maxzs[i];}
    maxz /= M;

    int Nx;
    double Lx, Ly, Lz, lenvoxel;
    vector<double> x_centres = vector<double>(Nvoxels);
    vector<double> y_centres_temp = vector<double>(Nvoxels);
    vector<double> z_centres_temp = vector<double>(Nvoxels);

    Nx = Nvoxels;
    Lx = maxx-minx;
    Ly = maxy-miny;
    Lz = maxz-minz;
    lenvoxel = Lx/Nvoxels;


    // Want to do something like this:
    //x_centres = 0.5*(2*ns+1)*lenvoxel+minx;
    for(int i=0; i<Nvoxels; i++){
        x_centres[i] = 0.5*(2.0*ns[i]+1.0)*lenvoxel+minx;
    }

    int centrevalue, edgevalue;
    int Ny      = 0;
    int Nz      = 0;
    int edgehit = 0;
    int counter = 0;
    double i     = 0;
    while (edgehit==0){
        centrevalue = 0.5*(2.0*i+1.0)*lenvoxel+miny;
        edgevalue   = lenvoxel*(i+1.0)+miny;
        y_centres_temp[counter] = centrevalue;
        i += 1.0;
        counter++;
        if ((maxy-centrevalue)<0.5*lenvoxel or abs(maxy-edgevalue)<0.1*lenvoxel){
            Ny = counter;
            edgehit = 1;
        }
    }
    i       = 0.0;
    counter = 0;
    edgehit = 0;
    while (edgehit==0){
        centrevalue = 0.5*(2.0*i+1.0)*lenvoxel+minz;
        edgevalue   = lenvoxel*(i+1.0)+minz;
        z_centres_temp[counter] = centrevalue;
        i += 1.0;
        counter++;
        if ((maxz-centrevalue)<0.5*lenvoxel or abs(maxz-edgevalue)<0.1*lenvoxel){
            Nz = counter;
            edgehit = 1;
        }
    }


    /* // I don't think I really need this. I only need to take care of Ny and Nz
    vector<double> y_centres = vector<double>(Ny);
    vector<double> z_centres = vector<double>(Nz);

    //y_centres = y_centres_temp[0:Ny];
    //z_centres = z_centres_temp[0:Nz];

    for(int i=0; i<Ny; i++){y_centres[i]=y_centres_temp[i]}
    */

    /* // This is the Python code for the base implementation of the voxellation
      // I guess I should put this in a function?
        def voxellate(Nx, Ny, Nz, x_centres, y_centres, z_centres, posvecs):
            voxmat   = np.zeros((Nx,Ny,Nz))
            n = 0
            for i in range(Nx):
                for j in range(Ny):
                    for k in range(Nz):
                        smallest_dist = 1e20 # No distance is this big.
                        xc            = x_centres[i]
                        yc            = y_centres[j]
                        zc            = z_centres[k]
                        centrevec     = np.array([xc,yc,zc])
                        # Loop over all the atoms.
                        for m in range(Nall):
                            vecthis = posvecs[m,:]
                            distvec = centrevec-vecthis       # This is broadcasting, and Numba likes that
                            dotprod = np.dot(distvec,distvec) # Numba likes numpy functions
                            if dotprod<smallest_dist:
                                 smallest_dist = dotprod
                        voxval         = np.sqrt(smallest_dist) # Another Numpy function
                        voxmat[i,j,k]  = voxval
                        n+=1
            return voxmat
    */

    // Put this in a function?
    //def voxellate(Nx, Ny, Nz, x_centres, y_centres, z_centres, posvecs):
    vector<vector<vector<double>>> voxelmat(Nx,vector<vector<double>>(Ny,vector <double>Nz));
    //vector<vector<vector<double>>> voxelmat; // = vector<vector<vector<double>>>(Nx,vector<vector<double>>(Ny,vector <double>Nz));
    int frame = 0;
    voxelmat = voxellate_basic(frame, Nx, Ny, Nz, x_centres, y_centres, z_centres, xes, ys, zs);



    /*
    */

    // ------------------- Voxellate ------------------- //


    cout << "Hello" << ". " << "It is me" << "\n" << totframes << "\n" << indices[0] << "\n" << indices[1] << endl;

    /*// A formulation like this would be handier
    char filename2[900];
    filename =
    sprintf(filename2, 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%.13f_Kbond%.12f_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj', M, N, gridspacing, wall, Kangle, Kbond, charge, T);
    */

    cout << filename << endl;

    return 0;
}

vector<int> int_linspace(int framefirst, int framelast, int Nframes) // This does not work for some reason...
{
    cout << "In int_linspace" << endl;
    vector<int> linspace_array = vector<int>(Nframes);
    int spacing = (framelast-framefirst)/(Nframes-1); // Double check
    for (int i; i<Nframes-1; i++)
    {
        linspace_array[i] = framefirst + i*spacing;
        cout << "i=" << i << "; value =" << framefirst + i*spacing << endl;
    }
    return linspace_array;
}

vector<vector<vector<double>>> voxellate_basic(int frame, int Nx, int Ny, int Nz, vector<double> x_centres, vector<double> y_centres, vector<double> z_centres, vector<vector<double>> xes, vector<vector<double>> ys, vector<vector<double>> zs)
{
    /* // This is the Python code for the base implementation of the voxellation
      // I guess I should put this in a function?
        def voxellate(Nx, Ny, Nz, x_centres, y_centres, z_centres, posvecs):
            voxmat   = np.zeros((Nx,Ny,Nz))
            n = 0
            for i in range(Nx):
                for j in range(Ny):
                    for k in range(Nz):
                        smallest_dist = 1e20 # No distance is this big.
                        xc            = x_centres[i]
                        yc            = y_centres[j]
                        zc            = z_centres[k]
                        centrevec     = np.array([xc,yc,zc])
                        # Loop over all the atoms.
                        for m in range(Nall):
                            vecthis = posvecs[m,:]
                            distvec = centrevec-vecthis       # This is broadcasting, and Numba likes that
                            dotprod = np.dot(distvec,distvec) # Numba likes numpy functions
                            if dotprod<smallest_dist:
                                 smallest_dist = dotprod
                        voxval         = np.sqrt(smallest_dist) # Another Numpy function
                        voxmat[i,j,k]  = voxval
                        n+=1
            return voxmat
    */

    double xc, yc, zc, xthis, ythis, zthis, smallest_dist;
    vector<double> diffvec       = vector<double>(3);
    vector<double> thisbeadvec   = vector<double>(3);
    vector<double> thiscentrevec = vector<double>(3);
    vector<vector<vector<double>>> voxmat(Nx,vector<vector<double>>(Ny,vector <double>Nz));
    // I'm gonna loop over the beads instead (there are fewer of them)
    for(int i=0; i<Nx; i++){
        xc = x_centres[i]; // Possibly use thiscentrevec directly. See what works best.
        thiscentrevec[0] = xc;
        for(int j=0; j<Ny; j++){
            yc = y_centres[i];
            thiscentrevec[1] = yc;
            for(int k=0; k<Nz; k++){
                zc = z_centres[i];
                thiscentrevec[2] = xc;
                smallest_dist = 1e20; // No distance is this big
                for(int m=0; m<Nall; m++){
                    xthis = xes[frame][m];
                    ythis = ys[frame][m];
                    zthis = zs[frame][m];
                    thisbeadvec[0] = xthis;
                    thisbeadvec[1] = ythis;
                    thisbeadvec[2] = zthis;
                    diffvec = thisbeadvec-thiscentrevec; // Do I really need to involve vectors, though? Why not simply use dx, dy, dz?

                }
            } // End loop over z-dir

        } // End loop over y-dir
    } // End loop over x-dir


}

