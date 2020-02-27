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
vector<vector<vector<double>>> voxellate_basic(int frame, int Nall, int Nx, int Ny, int Nz, vector<double> x_centres, vector<double> y_centres, vector<double> z_centres, vector<vector<double>> xes, vector<vector<double>> ys, vector<vector<double>> zs); // voxellate(Nx, Ny, Nz, x_centres, y_centres, z_centres, posvecs)

//vector<vector<double>> voxellate_basic(int Nx, int Ny, int Nz, vector<double> x_centres, vector<double> y_centres, vector<double> z_centres, vector<vector<double>> xes, vector<vector<double>> ys, vector<vector<double>> zs); // voxellate(Nx, Ny, Nz, x_centres, y_centres, z_centres, posvecs)

int main()
{
    int M, N, gridspacing;
    double wall, charge, T, Kangle, Kbond;
    string infilenamePrefix, outfilenamePrefix; //filename, outfilename, xfilename, yfilename, zfilename, datafilename;

    int Nvoxels    = 55;
    int totframes  = 10000000/10000;

    // File selection
    //string
    //filename = "chaingrid_quadratic_M9N101_gridspacing40_Langevin_wall1.042_Kangle14.0186574854529_Kbond140.186574854529_debye_kappa1_debyecutoff3_charge-1_T3_theta0is180_twofirst_are_fixed.lammpstrj";
    //filename = "chaingrid_quadratic_M9N101_gridspacing40_Langevin_wall1.042_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_T3_theta0is180_twofirst_are_fixed.lammpstrj";
    //filename = "C:/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/chaingrid_quadratic_M9N101_gridspacing40_Langevin_wall1.042_Kangle20_Kbond200_debye_kappa1_debyecutoff3_charge-1_T3_theta0is180_twofirst_are_fixed.lammpstrj";

    // File name for testing (look more into which files I've deleted and which I haven't later): # And look into paths.
    infilenamePrefix  = "/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/chaingrid_quadratic_M9N101_ljdebye1.042_angle_Langevin_wall1.042_Kangle25_Kbond2000_T310_theta0is180_cutoffs3_sigma1_firstatomfixed";
    outfilenamePrefix = "/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Fast_voxellation/chaingrid_quadratic_M9N101_ljdebye1.042_angle_Langevin_wall1.042_Kangle25_Kbond2000_T310_theta0is180_cutoffs3_sigma1_firstatomfixed";
    /*
    filename = "/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/chaingrid_quadratic_M9N101_ljdebye1.042_angle_Langevin_wall1.042_Kangle25_Kbond2000_T310_theta0is180_cutoffs3_sigma1_firstatomfixed.lammpstrj";
    outfilename = "/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Fast_voxellation/chaingrid_quadratic_M9N101_ljdebye1.042_angle_Langevin_wall1.042_Kangle25_Kbond2000_T310_theta0is180_cutoffs3_sigma1_firstatomfixed_step%i.txt";
    xfilename   = "/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Fast_voxellation/x_centres.txt";   // Strictly speaking, I don't need these, but they do make things a bit easier...
    yfilename   = "/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Fast_voxellation/y_centres.txt";
    zfilename   = "/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Fast_voxellation/z_centres.txt";
    datafilename = "/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Fast_voxellation/data.txt";       // This is basically all the information I need to reconstruct the data // But easier if I write bincenters to file.
    */

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
    char *filename = new char[100000]; // Possibly a long file name
    sprintf(filename, "%s.lammpstrj", infilenamePrefix.c_str());
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
    vector<double> y_centres = vector<double>(Nvoxels);
    vector<double> z_centres = vector<double>(Nvoxels);

    Nx = Nvoxels;
    Lx = maxx-minx;
    Ly = maxy-miny;
    Lz = maxz-minz;
    lenvoxel = Lx/Nvoxels;


    // Want to do something like this:
    //x_centres = 0.5*(2*ns+1)*lenvoxel+minx;
    for(int i=0; i<Nvoxels; i++){
        x_centres[i] = 0.5*(2.0*ns[i]+1.0)*lenvoxel+minx;
        //cout << "x_centres[" << i << "]:" << x_centres[i] << endl;
    }

    float centrevalue, edgevalue;
    int Ny      = 0;
    int Nz      = 0;
    int edgehit = 0;
    int counter = 0;
    double i    = 0;
    while (edgehit==0){
        centrevalue = 0.5*(2.0*i+1.0)*lenvoxel+miny;
        edgevalue   = lenvoxel*(i+1.0)+miny;
        y_centres[counter] = centrevalue;
        //cout << "y_centres[" << counter << "]:" << y_centres[counter] << "= " << centrevalue << endl;
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
        z_centres[counter] = centrevalue;
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


    for(int i=0; i<Ny; i++){y_centres[i]=y_centres[i]}
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

    // ------------------- Voxellate ------------------- //

    // Put this in a function?
    //def voxellate(Nx, Ny, Nz, x_centres, y_centres, z_centres, posvecs):
    vector<vector<vector<double>>> voxelmat = vector<vector<vector<double> > > (Nx,vector<vector<double>>(Ny,vector <double>(Nz)));
    //vector<vector<vector<double>>> voxelmat(Nx,vector<vector<double>>(Ny,vector <double>Nz));
    // vector<vector<vector<int> > > v = vector<vector<vector<int> > >( 5, vector<vector<int> >(3, vector<int>(2, 4)))
    //vector<vector<vector<double>>> voxelmat; // = vector<vector<vector<double>>>(Nx,vector<vector<double>>(Ny,vector <double>Nz));

    /*
    char *filenamez = new char[1000];                                // File name can have max 1000 characters
    sprintf(filenamez, "%s_spincorrelationfunctionz.txt", filenamePrefix.c_str() );   // Create filename with prefix and ending
    spcorFilez.open(filenamez);
    delete filenamez;
    */

    int frame;
    ofstream outFile;
    for(int l=0; l<Nout; l++){
        frame = outframes[l];
        cout << "l = " << l << ", frame " << frame << endl;
        voxelmat = voxellate_basic(l, Nall, Nx, Ny, Nz, x_centres, y_centres, z_centres, xes, ys, zs);

        //cout << "Voxellation complete" << endl;

        // File name
        char *ofilename = new char[100000]; // Possibly a long file name
        sprintf(ofilename, "/home/kine/Projects_PhD/P2_PolymerMD/Planar_brush/Fast_voxellation/chaingrid_quadratic_M9N101_ljdebye1.042_angle_Langevin_wall1.042_Kangle25_Kbond2000_T310_theta0is180_cutoffs3_sigma1_firstatomfixed_step%i", frame);
        //ofilename = sprintf(ofilename, "%s.txt", outfilenamePrefix.c_str());
        outFile.open(ofilename);
        delete ofilename;
        // Write data to file

        //cout << "Writing voxelmat to file" << endl;
        for(int i=0; i<Nx; i++){
            //cout << "l = " << l << ", frame: " << frame << ": i = " << i << " of " << Nx << endl;
            for(int j=0; j<Ny; j++){
                for(int k=0; k<Nz; k++){
                    // Write to file
                    outFile << std::setprecision(std::numeric_limits<double>::digits10+1) << voxelmat[i][j][k] << " ";
                } // Close loop over z
                outFile << endl;
            } // Close loop over y
        } // Close loop over x
        outFile.close();
        //cout << "Closing file" << endl;
        //cout << "outframes[0]: " << outframes[0] << ", outframes[1]: " << outframes[1] << ", outframes[2]: " << outframes[2] << endl;
    } // Close loop over frames


    //cout << "Done with writing voxelmat to file" << endl;

    ofstream dataFile;
    char *dfilename = new char[100000]; // Possibly a long file name
    sprintf(dfilename, "%s_data.txt", outfilenamePrefix.c_str());
    dataFile.open(dfilename); // Does this work?
    delete dfilename;
    dataFile << "len_voxel: " << lenvoxel << endl;
    dataFile << "Nx: " << Nx << " Ny: " << Ny << " Nz: " << Nz << endl;
    dataFile << "maxx: " << maxx << " minx: " << minx << endl;
    dataFile << "maxy: " << maxy << " miny: " << miny << endl;
    dataFile << "maxz: " << maxz << " minz: " << minz << endl;
    dataFile.close();

    // Voxel centre coordinates
    ofstream xFile;
    char *xfilename = new char[100000]; // Possibly a long file name
    sprintf(xfilename, "%s_x_centres.txt", outfilenamePrefix.c_str());
    xFile.open(xfilename);
    delete xfilename;

    ofstream yFile;
    char *yfilename = new char[100000]; // Possibly a long file name
    sprintf(yfilename, "%s_y_centres.txt", outfilenamePrefix.c_str());
    yFile.open(yfilename);
    delete yfilename;

    ofstream zFile;
    char *zfilename = new char[100000]; // Possibly a long file name
    sprintf(zfilename, "%s_z_centres.txt", outfilenamePrefix.c_str());
    zFile.open(zfilename);
    delete zfilename;

    // Voxel indices
    ofstream xindFile;
    char *xindfilename = new char[100000]; // Possibly a long file name
    sprintf(xindfilename, "%s_xind.txt", outfilenamePrefix.c_str());
    xindFile.open(xindfilename);
    delete xindfilename;

    ofstream yindFile;
    char *yindfilename = new char[100000]; // Possibly a long file name
    sprintf(yindfilename, "%s_yind.txt", outfilenamePrefix.c_str());
    yindFile.open(yindfilename);
    delete yindfilename;

    ofstream zindFile;
    char *zindfilename = new char[100000]; // Possibly a long file name
    sprintf(zindfilename, "%s_zind.txt", outfilenamePrefix.c_str());
    zindFile.open(zindfilename);
    delete zindfilename;

    for(int i=0; i<Nx; i++){
        for(int j=0; j<Ny; j++){
            for(int k=0; k<Nz; k++){
                // Write to file
                // Coordinates
                xFile << std::setprecision(std::numeric_limits<double>::digits10+1) << x_centres[i] << " "; // Is this really necessary?
                yFile << std::setprecision(std::numeric_limits<double>::digits10+1) << y_centres[j] << " ";
                zFile << std::setprecision(std::numeric_limits<double>::digits10+1) << z_centres[k] << " ";
                // Indices
                xindFile << i << " "; // Is this really necessary?
                yindFile << j << " ";
                zindFile << k << " ";
                //if(i==1 && j==1 && k==1){
                //    cout << "x_centres[i]:" <<  x_centres[i] << "; y_centres[j]:" <<  y_centres[i] << "; z_centres[k]:" <<  z_centres[i] << endl;
                //}
            } // Close loop over z
            // End lines, coordinate files
            xFile << endl;
            yFile << endl;
            zFile << endl;
            // End lines, index files
            xindFile << endl;
            yindFile << endl;
            zindFile << endl;
        } // Close loop over y
    } // Close loop over x

    //cout << "Should have written to coordinates files now" << endl;
    // Close coordinate files
    xFile.close();
    yFile.close();
    zFile.close();
    // Close index files
    xindFile.close();
    yindFile.close();
    zindFile.close();

    /*
    */

    /*// A formulation like this would be handier
    char filename2[900];
    sprintf(filename2, 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%.13f_Kbond%.12f_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj', M, N, gridspacing, wall, Kangle, Kbond, charge, T);
    */

    //cout << "outframes[0]: " << outframes[0] << ", outframes[1]: " << outframes[1] << ", outframes[2]: " << outframes[2] << endl;

    for(int i=0; i<Nout; i++){cout << "outframes[" << i << "]:" << outframes[i] << " ";}
    cout << endl;

    cout << "File done" << endl;

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

vector<vector<vector<double>>> voxellate_basic(int frame, int Nall, int Nx, int Ny, int Nz, vector<double> x_centres, vector<double> y_centres, vector<double> z_centres, vector<vector<double>> xes, vector<vector<double>> ys, vector<vector<double>> zs)
{
    //cout << "In voxellate_basic" << endl;
    double dx, dy, dz, xc, yc, zc, xthis, ythis, zthis, smallest_dist, dotprod;
    //cout << "Variables declared" << endl;
    vector<vector<vector<double>>> voxmat = vector<vector<vector<double>>> (Nx,vector<vector<double>>(Ny,vector <double>(Nz)));
    //cout << "Matrix initiated" << endl;
    // vector<vector<vector<int> > > v = vector<vector<vector<int> > >( 5, vector<vector<int> >(3, vector<int>(2, 4)))
    // I'm gonna loop over the beads instead (there are fewer of them)
    for(int i=0; i<Nx; i++){
        //cout << "In loop in voxellate_basic, i = " << i << endl;
        xc = x_centres[i];
        //cout << "x_centres[" << i << "] extracted" << endl;
        for(int j=0; j<Ny; j++){
            //cout << "In loop in voxellate_basic, i = " << i << ", j = " << j << endl;
            yc = y_centres[j];
            //cout << "y_centres[" << j << "] extracted" << endl;
            for(int k=0; k<Nz; k++){
                //cout << "In loop in voxellate_basic, i = " << i << ", j = " << j << ", k = " << k << endl;
                zc = z_centres[k];
                //cout << "z_centres[" << k << "] extracted" << endl;
                smallest_dist = 1e20; // No distance is this big
                for(int m=0; m<Nall; m++){
                    //cout << "In loop in voxellate_basic, i = " << i << ", j = " << j << ", k = " << k << ", m = " << m << endl;
                    xthis = xes[frame][m];
                    ythis = ys[frame][m];
                    zthis = zs[frame][m];
                    //cout << "xthis etc. extracted" << endl;
                    dx = xthis-xc;
                    dy = ythis-yc;
                    dz = zthis-zc;
                    //cout << "dx, etc. found" << endl;
                    dotprod = dx*dx +dy*dy +dz*dz;
                    //cout << "dotprod calculated" << endl;
                    if(dotprod<smallest_dist){smallest_dist=dotprod;}
                    //cout << "dotprod tested for smalles distance" << endl;
                }
                //cout << "Setting voxmat element" << endl;
                voxmat[i][j][k] = sqrt(smallest_dist);
                //cout << "voxmat element set" << endl;
            } // End loop over z-dir
        } // End loop over y-dir
    } // End loop over x-dir

    //cout << "Loop done. Returning in a moment" << endl;
    return voxmat;
}

