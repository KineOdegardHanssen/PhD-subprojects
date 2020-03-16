#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <cmath>
#include <string>
//#include <cstdio>
#include<stdio.h>
#include<sstream>

using namespace std;

// The preferred functions, I'd say: // Minor adjustment: intsigma gives the radius/half length of the structure, i.e. intsigma = 2 gives a sphere of diameter 5. Makes it easier to locate midpoint
// Non-periodic BCs
void randomwalk(int sigma, int d, int Nblocks, int Nrun, int Nrealizations, int Nsections);
void matrixwalk_hard_easyinput(int sigma, int d, int Nblocks, int Nrun, int Nrealizations, int Nsections, bool printmatrix);
void matrixmc_hard_easyinput(int sigma, int d, int Nblocks, int Nrun, int Nrealizations, int Nsections, bool printmatrix);
void matrixmc_potential_easyinput(int intsigma, int intd, int Nblocks, int Nrun, int Nrealizations, int Nsections, double beta, double sigma, double power, double soften); // Power is expensive
// Periodic BCs
void randomwalk_PBC(int sigma, int d, int Nblocks, int Nrun, int Nrealizations, int Nsections, int printevery, bool pbc_codetesting); // Could have made this totally open too... (The current form made sense when I didn't have PBCs.)
void matrixwalk_hard_easyinput_PBC(int sigma, int d, int Nblocks, int Nrun, int Nrealizations, int Nsections, int printevery, bool printmatrix, bool pbc_codetesting);
void matrixmc_hard_easyinput_PBC(int sigma, int d, int Nblocks, int Nrun, int Nrealizations, int Nsections, int printevery, bool printmatrix, bool pbc_codetesting);
void matrixmc_potential_easyinput_PBC(int intsigma, int intd, int Nblocks, int Nrun, int Nrealizations, int Nsections, int printevery, double beta, double sigma, double power, double soften, bool pbc_codetesting); // Power is expensive
void nearwall_rw(int maxdistance, int Nrun, int Nrealizations, int Nsections, int printevery); // No structure besides the wall. A totally open system
void nearwall_mc(int maxdistance, int Nrun, int Nrealizations, int Nsections, int printevery); // No structure besides the wall. A totally open system
void randomwalk_1D(int Nrun, int Nrealizations, int Nsections, int printevery);

// I could consider removing these functions:
void matrixwalk_hard(int Nblocks, int blocksize, int Nrun, int Nrealizations, int Nsections);
void matrixmc_hard(int Nblocks, int blocksize, int Nrun, int Nrealizations, int Nsections);
void matrixmc_potential(int Nblocks, int blocksize, int Nrun, int Nrealizations, int Nsections, double beta, double sigma, double power, double soften); // Power is expensive

// Auxilliary functions:
vector<int> give_startpoints(int Nrun, int Nsections);
double give_energy(int thisx, int thisy, int Npots, double sigma, double power, double soften, vector<double> x_centres, vector<double> y_centres);

int main()
{
    int intd, intsigma, Nrun, Nblocks, Nrealizations, Nsections, blocksize, maxdistance, printevery;//Lx, Ly;
    double beta, sigma, power, soften;
    bool printmatrix, pbc_codetesting;

    // All
    Nrun = 20000;
    Nblocks = 3;
    Nrealizations = 1000;
    printevery    = 100;  // Print every this many time steps

    // Blocky implementation (spacing and obstacles are the same length)
    blocksize = 2;
    Nsections = 5;

    // Flexible hard-block implementation (obstacles and spacing can be of different lengths)
    intd = 1;
    intsigma = 10;
    printmatrix = true;

    // Potential
    beta   = 3.0;
    sigma  = 10;
    power  = 6;
    soften = 1;

    // For walk/mc near wall
    maxdistance = 10;

    // PBC
    pbc_codetesting = false;

    //----- Function calls -----//
    //vector<int> startpoints = give_startpoints(Nrun, Nsections); // Call this inside function instead
    //matrixwalk_hard(Nblocks, blocksize, Nrun, Nrealizations, Nsections);
    //matrixmc_hard(Nblocks, blocksize, Nrun, Nrealizations, Nsections);
    //matrixmc_potential(Nblocks, blocksize, Nrun, Nrealizations, Nsections, beta, sigma, power, soften);
    //matrixwalk_hard_easyinput(intsigma, intd, Nblocks, Nrun, Nrealizations, Nsections, printmatrix);
    //matrixmc_potential_easyinput(intsigma, intd, Nblocks, Nrun, Nrealizations, Nsections, beta, sigma, power, soften); // Is there any
    //matrixmc_hard_easyinput(intsigma, intd, Nblocks, Nrun, Nrealizations, Nsections, printmatrix);
    //randomwalk(intsigma, intd, Nblocks, Nrun, Nrealizations, Nsections);
    //-- PBC --//
    //randomwalk_PBC(intsigma, intd, Nblocks, Nrun, Nrealizations, Nsections, printevery, pbc_codetesting);
    //matrixwalk_hard_easyinput_PBC(intsigma, intd, Nblocks, Nrun, Nrealizations, Nsections, printevery, printmatrix, pbc_codetesting);
    //matrixmc_hard_easyinput_PBC(intsigma, intd, Nblocks, Nrun, Nrealizations, Nsections, printevery, printmatrix, pbc_codetesting);
    //matrixmc_potential_easyinput_PBC(intsigma, intd, Nblocks, Nrun, Nrealizations, Nsections, printevery, beta, sigma, power, soften, pbc_codetesting); // Power is expensive
    //nearwall_rw(maxdistance, Nrun, Nrealizations, Nsections, printevery);
    //nearwall_mc(maxdistance, Nrun, Nrealizations, Nsections, printevery);
    randomwalk_1D(Nrun, Nrealizations, Nsections, printevery);

    cout << "Hello World!" << endl;
    return 0;
}

//---------------- Flexible functions ----------------//
//-- Non-periodic --//
void randomwalk(int sigma, int d, int Nblocks, int Nrun, int Nrealizations, int Nsections)
{   // Assuming quadratic matrix
    // Sigma: 'radius' of obstactle (i.e. half the length since it is quadratic)
    //---- Matrix set up ----//
    int Nmat, blocksize;
    blocksize = 2*sigma+1;
    Nmat = (Nblocks+1)*d + Nblocks*blocksize; // Size of matrix: contribution from spacings + contributions from obstacles (Structure: sp-obst-sp-obst-...-sp-obst)

    //---- Partitioning set up ----//
    vector<int> startpoints = give_startpoints(Nrun, Nsections);
    int len_sections = Nrun/Nsections; // I have already done this. Maybe having the function is not that handy...

    //---- Set up of run ----//
    // Random number distributions
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,Nmat-1);
    std::uniform_int_distribution<int> steps_distr(0,1); // Move: +/-1. Should I be able to keep still too?
    std::uniform_int_distribution<int> dir_distr(0,1);

    // Variables
    bool moveit;
    //double thisR2;    // Is this a conflict? Should this be an int?
    int starti, startj, indi, indj, nextindi, nextindj, dir, step, thisR2, nx, ny;
    vector<double> walk_R2(Nrun);                                           // Average R^2
    vector<double> walk_R2_rms(Nrun);                                       // RMS R^2
    vector<vector<double>> walk_R2_store(Nrealizations, vector<double>(Nrun)); // For calculation of the standard deviation
    vector<vector<double>> walk_x(Nrealizations, vector<double>(Nrun));        // Useful in case I want to look at the walks for different starting points (post-processing)
    vector<vector<double>> walk_y(Nrealizations, vector<double>(Nrun));
    vector<vector<double>> walk_R2_sections(Nsections, vector<double>(len_sections));                                                      // Sections
    vector<vector<double>> walk_R2_sections_rms(Nsections, vector<double>(len_sections));                                                  // Sections, rms
    vector<vector<vector<double>>> walk_R2_sections_store(Nrealizations, vector<vector<double>>(Nsections, vector<double>(len_sections))); // For calculation of the standard deviation


    //---- Performing walk ----//
    nx = 0; ny = 0;
    for(int i=0; i<Nrealizations; i++){
        // Start the walk by placing the bead in a random position
        // and check that the position is not occupied
        indi = distribution(generator); //location();//indi = rand() % Nmat; // Random numbers this way (from <random>) is a bit better than rand() from cstdlib, it seems.
        indj = distribution(generator); //location();//indj = rand() % Nmat;
        starti = indi;
        startj = indj;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl;
        for(int j=0; j<Nrun; j++){
            moveit = false;
            while(!moveit){    // Should I let the walker stand still if the move is not accepted?
                step = 2*steps_distr(generator)-1;
                //cout << "step: " << step << endl;
                dir = dir_distr(generator);
                //cout << "dir: " << dir << endl;
                if(dir==0){
                    //cout << "dir x chosen" << endl;
                    nextindi = indi + step;
                    nextindj = indj;
                    if(nextindi<Nmat && nextindi>-1){
                        moveit=true;
                        nx++;
                    }
                }
                else{
                    //cout << "dir y chosen" << endl;
                    nextindi = indi;
                    nextindj = indj + step;
                    if(nextindj<Nmat && nextindj>-1){ // If the site exists
                        moveit=true;
                        ny++;
                    }
                }
            }
            // Move the walker
            indi = nextindi;
            indj = nextindj;
            // store data somewhere
            // Average moves
            thisR2 = (indi-starti)*(indi-starti) + (indj-startj)*(indj-startj);
            walk_R2[j] += thisR2;
            walk_R2_store[i][j] = thisR2;
            // Saving coordinates
            walk_x[i][j] = indi;
            walk_y[i][j] = indj;
        }
        //cout << "indi: " << indi << ", indj: " << indj << endl;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl << "-------------------------------------" << endl;
    }

    for(int i=0; i<Nrun; i++){
        walk_R2[i]/=Nrealizations;
        //cout << "walk_R2[" << i << "]: " << walk_R2[i] << endl;
    } // Averaging complete


    // Finding rms values
    for(int j=0; j<Nrun; j++){
        for(int i=0; i<Nrealizations; i++){
            walk_R2_rms[j] += (walk_R2_store[i][j]-walk_R2[j])*(walk_R2_store[i][j]-walk_R2[j]);
        }
        walk_R2_rms[j] = sqrt(walk_R2_rms[j]/(Nrealizations-1));
    }

    //cout << "nx: " << nx << endl << "ny: " << ny << endl;
    //cout << "Nmat-1: " << Nmat-1 << endl;


    //----- Partitioning -----//

    int startpoint;
    double startx, starty, currx, curry, dx, dy, dotprod;
    for(int i=0; i<Nrealizations; i++){
        for(int j=0; j<Nsections; j++){
            // Extract starting points
            startpoint = startpoints[j];
            startx = walk_x[i][startpoint];
            starty = walk_y[i][startpoint];
            for(int k=1; k<len_sections; k++){ // Start at 1 since first point in array will be 0 anyways
                currx   = walk_x[i][startpoint+k];
                curry   = walk_y[i][startpoint+k];
                dx      = currx-startx;
                dy      = curry-starty;
                dotprod = dx*dx+dy*dy;
                walk_R2_sections[j][k]          += dotprod;
                walk_R2_sections_store[i][j][k]  = dotprod;
            } // End loop over steps (len_sections)
        } // End loop over sections
    } // End loop over realizations

   // Averaging
    for(int i = 0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            walk_R2_sections[i][j] /= Nrealizations;
        } // End loop over steps (len_sections)
    } // End loop over sections

    // rms
    for(int i=0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            for(int k=0; k<Nrealizations; k++){
                walk_R2_sections_rms[i][j] = (walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j])*(walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j]);
            } // End loop over realizations
            walk_R2_sections_rms[i][j] = sqrt(walk_R2_sections_rms[i][j]/(Nrealizations-1));
        } // End loop over steps (len_sections)
    } // End loop over sections


    //----- Writing to file -----//
    //string outfilenamePrefix;


    ofstream outFile;
    char *filename = new char[100000]; // Possibly a long file name
    sprintf(filename, "Nsteps%i_Nreal%i_sigma%i_d%i_Nblocks%i_randomwalk_R2_basic_mc", Nrun, Nrealizations, sigma, d, Nblocks);
    outFile.open(filename);
    delete filename;


    for(int i=0; i<Nrun; i++){
        outFile << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl;
        //outFile << std::setprecision(std::numeric_limits<double>::digits10+1) << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl; // Is this long printing really neccessary? // Possibly (probably?) gonna use this later.
    }
    outFile.close();

    ofstream poutFile;
    char *pfilename = new char[100000]; // Possibly a long file name
    sprintf(pfilename, "Nsteps%i_Nreal%i_sigma%i_d%i_Nblocks%i_Npart%i_randomwalk_R2_basic_mc", Nrun, Nrealizations, sigma, d, Nblocks, Nsections);
    poutFile.open(pfilename);
    delete pfilename;

    for(int i=0; i<Nsections; i++){
        poutFile << "--- Section " << i <<" ---" << endl;
        for(int j=0; j<len_sections; j++){
            poutFile << j << " " << walk_R2_sections[i][j] << " " << walk_R2_sections_rms[i][j] << endl;
        } // End loop over steps (len_sections)
    } // End loop over sections
    poutFile.close();
}



void matrixwalk_hard_easyinput(int sigma, int d, int Nblocks, int Nrun, int Nrealizations, int Nsections, bool printmatrix)
{   // Assuming quadratic matrix
    // Sigma: 'radius' of obstactle (i.e. half the length since it is quadratic)
    //---- Matrix set up ----//
    int Nmat, blocksize;
    blocksize = 2*sigma+1;
    Nmat = (Nblocks+1)*d + Nblocks*blocksize; // Size of matrix: contribution from spacings + contributions from obstacles (Structure: sp-obst-sp-obst-...-sp-obst)
    vector<vector<int>> walkmat(Nmat, vector<int>(Nmat));

    // Setting matrix blocks
    // Must be a handier way to do this?:
    for(int i=0; i<Nblocks; i++){
        for(int j=0; j<Nblocks; j++){
            for(int k=0; k<blocksize; k++){
                for(int l=0; l<blocksize; l++){
                    walkmat[(i+1)*d+i*blocksize+k][(j+1)*d+j*blocksize+l] = 1;
                }
            }
        }
    }

    //---- Partitioning set up ----//
    vector<int> startpoints = give_startpoints(Nrun, Nsections);
    int len_sections = Nrun/Nsections; // I have already done this. Maybe having the function is not that handy...

    //---- Set up of run ----//
    // Random number distributions
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,Nmat-1);
    std::uniform_int_distribution<int> steps_distr(0,1); // Move: +/-1. Should I be able to keep still too?
    std::uniform_int_distribution<int> dir_distr(0,1);

    // Variables
    bool free, moveit;
    //double thisR2;    // Is this a conflict? Should this be an int?
    int starti, startj, indi, indj, nextindi, nextindj, dir, step, thisR2, nx, ny;
    vector<double> walk_R2(Nrun);                                           // Average R^2
    vector<double> walk_R2_rms(Nrun);                                       // RMS R^2
    vector<vector<double>> walk_R2_store(Nrealizations, vector<double>(Nrun)); // For calculation of the standard deviation
    vector<vector<double>> walk_x(Nrealizations, vector<double>(Nrun));        // Useful in case I want to look at the walks for different starting points (post-processing)
    vector<vector<double>> walk_y(Nrealizations, vector<double>(Nrun));
    vector<vector<double>> walk_R2_sections(Nsections, vector<double>(len_sections));                                                      // Sections
    vector<vector<double>> walk_R2_sections_rms(Nsections, vector<double>(len_sections));                                                  // Sections, rms
    vector<vector<vector<double>>> walk_R2_sections_store(Nrealizations, vector<vector<double>>(Nsections, vector<double>(len_sections))); // For calculation of the standard deviation


    //---- Performing walk ----//
    nx = 0; ny = 0;
    for(int i=0; i<Nrealizations; i++){
        // Start the walk by placing the bead in a random position
        // and check that the position is not occupied
        free = false;
        while(!free){
            indi = distribution(generator); //location();//indi = rand() % Nmat; // Random numbers this way (from <random>) is a bit better than rand() from cstdlib, it seems.
            indj = distribution(generator); //location();//indj = rand() % Nmat;
            //cout << "walkmat["<<indi<<"]["<<indj<< "]: " << walkmat[indi][indj] << endl;
            if(walkmat[indi][indj]==0){free=true;} // Only place the bead on a free site
        }
        starti = indi;
        startj = indj;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl;
        for(int j=0; j<Nrun; j++){
            moveit = false;
            while(!moveit){    // Should I let the walker stand still if the move is not accepted?
                step = 2*steps_distr(generator)-1;
                //cout << "step: " << step << endl;
                dir = dir_distr(generator);
                //cout << "dir: " << dir << endl;
                if(dir==0){
                    //cout << "dir x chosen" << endl;
                    nextindi = indi + step;
                    nextindj = indj;
                    if(nextindi<Nmat && nextindi>-1){
                        if(walkmat[nextindi][nextindj]==0){
                            moveit=true;
                            nx++;
                        }
                    }
                }
                else{
                    //cout << "dir y chosen" << endl;
                    nextindi = indi;
                    nextindj = indj + step;
                    if(nextindj<Nmat && nextindj>-1){ // If the site exists
                        if(walkmat[nextindi][nextindj]==0){
                            moveit=true;
                            ny++;
                        }
                    }
                }
            }
            // Move the walker
            indi = nextindi;
            indj = nextindj;
            // store data somewhere
            // Average moves
            thisR2 = (indi-starti)*(indi-starti) + (indj-startj)*(indj-startj);
            walk_R2[j] += thisR2;
            walk_R2_store[i][j] = thisR2;
            // Saving coordinates
            walk_x[i][j] = indi;
            walk_y[i][j] = indj;
        }
        //cout << "indi: " << indi << ", indj: " << indj << endl;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl << "-------------------------------------" << endl;
    }

    for(int i=0; i<Nrun; i++){
        walk_R2[i]/=Nrealizations;
        //cout << "walk_R2[" << i << "]: " << walk_R2[i] << endl;
    } // Averaging complete


    // Finding rms values
    for(int j=0; j<Nrun; j++){
        for(int i=0; i<Nrealizations; i++){
            walk_R2_rms[j] += (walk_R2_store[i][j]-walk_R2[j])*(walk_R2_store[i][j]-walk_R2[j]);
        }
        walk_R2_rms[j] = sqrt(walk_R2_rms[j]/(Nrealizations-1));
    }

    //cout << "nx: " << nx << endl << "ny: " << ny << endl;
    //cout << "Nmat-1: " << Nmat-1 << endl;


    //----- Partitioning -----//

    int startpoint;
    double startx, starty, currx, curry, dx, dy, dotprod;
    for(int i=0; i<Nrealizations; i++){
        for(int j=0; j<Nsections; j++){
            // Extract starting points
            startpoint = startpoints[j];
            startx = walk_x[i][startpoint];
            starty = walk_y[i][startpoint];
            for(int k=1; k<len_sections; k++){ // Start at 1 since first point in array will be 0 anyways
                currx   = walk_x[i][startpoint+k];
                curry   = walk_y[i][startpoint+k];
                dx      = currx-startx;
                dy      = curry-starty;
                dotprod = dx*dx+dy*dy;
                walk_R2_sections[j][k]          += dotprod;
                walk_R2_sections_store[i][j][k]  = dotprod;
            } // End loop over steps (len_sections)
        } // End loop over sections
    } // End loop over realizations

   // Averaging
    for(int i = 0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            walk_R2_sections[i][j] /= Nrealizations;
        } // End loop over steps (len_sections)
    } // End loop over sections

    // rms
    for(int i=0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            for(int k=0; k<Nrealizations; k++){
                walk_R2_sections_rms[i][j] = (walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j])*(walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j]);
            } // End loop over realizations
            walk_R2_sections_rms[i][j] = sqrt(walk_R2_sections_rms[i][j]/(Nrealizations-1));
        } // End loop over steps (len_sections)
    } // End loop over sections


    //----- Writing to file -----//
    //string outfilenamePrefix;


    ofstream outFile;
    char *filename = new char[100000]; // Possibly a long file name
    sprintf(filename, "Nsteps%i_Nreal%i_sigma%i_d%i_Nblocks%i_hardpot_R2_basic_mc", Nrun, Nrealizations, sigma, d, Nblocks);
    outFile.open(filename);
    delete filename;


    for(int i=0; i<Nrun; i++){
        outFile << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl;
        //outFile << std::setprecision(std::numeric_limits<double>::digits10+1) << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl; // Is this long printing really neccessary? // Possibly (probably?) gonna use this later.
    }
    outFile.close();

    ofstream poutFile;
    char *pfilename = new char[100000]; // Possibly a long file name
    sprintf(pfilename, "Nsteps%i_Nreal%i_sigma%i_d%i_Nblocks%i_Npart%i_hardpot_R2_basic_mc", Nrun, Nrealizations, sigma, d, Nblocks, Nsections);
    poutFile.open(pfilename);
    delete pfilename;

    for(int i=0; i<Nsections; i++){
        poutFile << "--- Section " << i <<" ---" << endl;
        for(int j=0; j<len_sections; j++){
            poutFile << j << " " << walk_R2_sections[i][j] << " " << walk_R2_sections_rms[i][j] << endl;
        } // End loop over steps (len_sections)
    } // End loop over sections
    poutFile.close();

    if(printmatrix){
        // Print matrix for testing:
        ofstream moutFile;
        char *mfilename = new char[100000]; // Possibly a long file name
        sprintf(mfilename, "sigma%i_d%i_Nblocks%i_Npart%i_matrix", Nrun, Nrealizations, sigma, d, Nblocks, Nsections);
        moutFile.open(mfilename);
        delete mfilename;

        for(int i=0; i<Nmat; i++){
            for(int j=0; j<Nmat; j++){
                moutFile << walkmat[i][j] << " ";
            } // End loop over steps (len_sections)
            moutFile << endl;
        } // End loop over sections
        moutFile.close();
    } // End if(printmatrix)

}

// MC
void matrixmc_hard_easyinput(int sigma, int d, int Nblocks, int Nrun, int Nrealizations, int Nsections,  bool printmatrix)
{   // Assuming quadratic matrix
    // Sigma: 'radius' of obstactle (i.e. half the length since it is quadratic)
    //---- Matrix set up ----//
    int Nmat, blocksize;
    blocksize = 2*sigma+1;
    Nmat = (Nblocks+1)*d + Nblocks*blocksize; // Size of matrix: contribution from spacings + contributions from obstacles (Structure: sp-obst-sp-obst-...-sp-obst)
    vector<vector<int>> walkmat(Nmat, vector<int>(Nmat));

    // Setting matrix blocks
    // Must be a handier way to do this?:
    for(int i=0; i<Nblocks; i++){
        for(int j=0; j<Nblocks; j++){
            for(int k=0; k<blocksize; k++){
                for(int l=0; l<blocksize; l++){
                    walkmat[(i+1)*d+i*blocksize+k][(j+1)*d+j*blocksize+l] = 1;
                }
            }
        }
    }

    //---- Partitioning set up ----//
    vector<int> startpoints = give_startpoints(Nrun, Nsections);
    int len_sections = Nrun/Nsections; // I have already done this. Maybe having the function is not that handy...

    //---- Set up of run ----//
    // Random number distributions
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,Nmat-1);
    std::uniform_int_distribution<int> steps_distr(0,1); // Move: +/-1. Should I be able to keep still too?
    std::uniform_int_distribution<int> dir_distr(0,1);

    // Variables
    bool free;
    //double thisR2;    // Is this a conflict? Should this be an int?
    int starti, startj, indi, indj, nextindi, nextindj, dir, step, thisR2, nx, ny;
    vector<double> walk_R2(Nrun);                                           // Average R^2
    vector<double> walk_R2_rms(Nrun);                                       // RMS R^2
    vector<vector<double>> walk_R2_store(Nrealizations, vector<double>(Nrun)); // For calculation of the standard deviation
    vector<vector<double>> walk_x(Nrealizations, vector<double>(Nrun));        // Useful in case I want to look at the walks for different starting points (post-processing)
    vector<vector<double>> walk_y(Nrealizations, vector<double>(Nrun));
    vector<vector<double>> walk_R2_sections(Nsections, vector<double>(len_sections));                                                      // Sections
    vector<vector<double>> walk_R2_sections_rms(Nsections, vector<double>(len_sections));                                                  // Sections, rms
    vector<vector<vector<double>>> walk_R2_sections_store(Nrealizations, vector<vector<double>>(Nsections, vector<double>(len_sections))); // For calculation of the standard deviation


    //---- Performing walk ----//
    nx = 0; ny = 0;
    for(int i=0; i<Nrealizations; i++){
        // Start the walk by placing the bead in a random position
        // and check that the position is not occupied
        free = false;
        while(!free){
            indi = distribution(generator); //location();//indi = rand() % Nmat; // Random numbers this way (from <random>) is a bit better than rand() from cstdlib, it seems.
            indj = distribution(generator); //location();//indj = rand() % Nmat;
            //cout << "walkmat["<<indi<<"]["<<indj<< "]: " << walkmat[indi][indj] << endl;
            if(walkmat[indi][indj]==0){free=true;} // Only place the bead on a free site
        }
        starti = indi;
        startj = indj;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl;
        for(int j=0; j<Nrun; j++){
            step = 2*steps_distr(generator)-1;
            //cout << "step: " << step << endl;
            dir = dir_distr(generator);
            //cout << "dir: " << dir << endl;
            if(dir==0){
                //cout << "dir x chosen" << endl;
                nextindi = indi + step;
                nextindj = indj;
                if(nextindi<Nmat && nextindi>-1){
                    if(walkmat[nextindi][nextindj]==0){
                        nx++;
                    }
                    else{
                        nextindi = indi;
                        nextindj = indj;
                    }
                }
                else{
                    nextindi = indi;
                    nextindj = indj;
                }
            }
            else{
                //cout << "dir y chosen" << endl;
                nextindi = indi;
                nextindj = indj + step;
                if(nextindj<Nmat && nextindj>-1){ // If the site exists
                    if(walkmat[nextindi][nextindj]==0){
                        ny++;
                    }
                    else{
                        nextindi = indi;
                        nextindj = indj;
                    }
                }
                else{
                    nextindi = indi;
                    nextindj = indj;
                }
            }
            // Move the walker
            indi = nextindi;
            indj = nextindj;
            // store data somewhere
            // Average moves
            thisR2 = (indi-starti)*(indi-starti) + (indj-startj)*(indj-startj);
            walk_R2[j] += thisR2;
            walk_R2_store[i][j] = thisR2;
            // Saving coordinates
            walk_x[i][j] = indi;
            walk_y[i][j] = indj;
        }
        //cout << "indi: " << indi << ", indj: " << indj << endl;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl << "-------------------------------------" << endl;
    }

    for(int i=0; i<Nrun; i++){
        walk_R2[i]/=Nrealizations;
        //cout << "walk_R2[" << i << "]: " << walk_R2[i] << endl;
    } // Averaging complete


    // Finding rms values
    for(int j=0; j<Nrun; j++){
        for(int i=0; i<Nrealizations; i++){
            walk_R2_rms[j] += (walk_R2_store[i][j]-walk_R2[j])*(walk_R2_store[i][j]-walk_R2[j]);
        }
        walk_R2_rms[j] = sqrt(walk_R2_rms[j]/(Nrealizations-1));
    }

    //cout << "nx: " << nx << endl << "ny: " << ny << endl;
    //cout << "Nmat-1: " << Nmat-1 << endl;


    //----- Partitioning -----//

    int startpoint;
    double startx, starty, currx, curry, dx, dy, dotprod;
    for(int i=0; i<Nrealizations; i++){
        for(int j=0; j<Nsections; j++){
            // Extract starting points
            startpoint = startpoints[j];
            startx = walk_x[i][startpoint];
            starty = walk_y[i][startpoint];
            for(int k=1; k<len_sections; k++){ // Start at 1 since first point in array will be 0 anyways
                currx   = walk_x[i][startpoint+k];
                curry   = walk_y[i][startpoint+k];
                dx      = currx-startx;
                dy      = curry-starty;
                dotprod = dx*dx+dy*dy;
                walk_R2_sections[j][k]          += dotprod;
                walk_R2_sections_store[i][j][k]  = dotprod;
            } // End loop over steps (len_sections)
        } // End loop over sections
    } // End loop over realizations

   // Averaging
    for(int i = 0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            walk_R2_sections[i][j] /= Nrealizations;
        } // End loop over steps (len_sections)
    } // End loop over sections

    // rms
    for(int i=0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            for(int k=0; k<Nrealizations; k++){
                walk_R2_sections_rms[i][j] = (walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j])*(walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j]);
            } // End loop over realizations
            walk_R2_sections_rms[i][j] = sqrt(walk_R2_sections_rms[i][j]/(Nrealizations-1));
        } // End loop over steps (len_sections)
    } // End loop over sections


    //----- Writing to file -----//
    //string outfilenamePrefix;


    ofstream outFile;
    char *filename = new char[100000]; // Possibly a long file name
    sprintf(filename, "Nsteps%i_Nreal%i_sigma%i_d%i_Nblocks%i_hardpot_R2_basic", Nrun, Nrealizations, sigma, d, Nblocks);
    outFile.open(filename);
    delete filename;


    for(int i=0; i<Nrun; i++){
        outFile << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl;
        //outFile << std::setprecision(std::numeric_limits<double>::digits10+1) << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl; // Is this long printing really neccessary? // Possibly (probably?) gonna use this later.
    }
    outFile.close();

    ofstream poutFile;
    char *pfilename = new char[100000]; // Possibly a long file name
    sprintf(pfilename, "Nsteps%i_Nreal%i_sigma%i_d%i_Nblocks%i_Npart%i_hardpot_R2_basic", Nrun, Nrealizations, sigma, d, Nblocks, Nsections);
    poutFile.open(pfilename);
    delete pfilename;

    for(int i=0; i<Nsections; i++){
        poutFile << "--- Section " << i <<" ---" << endl;
        for(int j=0; j<len_sections; j++){
            poutFile << j << " " << walk_R2_sections[i][j] << " " << walk_R2_sections_rms[i][j] << endl;
        } // End loop over steps (len_sections)
    } // End loop over sections
    poutFile.close();

    if(printmatrix){
        // Print matrix for testing:
        ofstream moutFile;
        char *mfilename = new char[100000]; // Possibly a long file name
        sprintf(mfilename, "sigma%i_d%i_Nblocks%i_Npart%i_matrix", Nrun, Nrealizations, sigma, d, Nblocks, Nsections);
        moutFile.open(mfilename);
        delete mfilename;

        for(int i=0; i<Nmat; i++){
            for(int j=0; j<Nmat; j++){
                moutFile << walkmat[i][j] << " ";
            } // End loop over steps (len_sections)
            moutFile << endl;
        } // End loop over sections
        moutFile.close();
    } // End if(printmatrix)

}


void matrixmc_potential_easyinput(int intsigma, int intd, int Nblocks, int Nrun, int Nrealizations, int Nsections, double beta, double sigma, double power, double soften)
{
    int Nmat, blocksize;
    blocksize = 2*intsigma+1;
    Nmat = (Nblocks+1)*intd + Nblocks*blocksize ; // Size of matrix: contribution from spacings + contributions from obstacles (Structure: sp-obst-sp-obst-...-sp-obst)
    //vector<vector<int>> walkmat(Nmat, vector<int>(Nmat)); // I don't really need walkmat, at least not yet.

    //---- Potential set up ----//
    // Number of 'obstacles' given by Nblocks x Nblocks
    // Centre them evenly in a grid-like manner
    int Npots = Nblocks*Nblocks;
    vector<double> centres_x = vector<double>(Npots);
    vector<double> centres_y = vector<double>(Npots);


    ofstream brushFile;
    char *bfilename = new char[100000]; // Possibly a long file name
    sprintf(bfilename, "intsigma%i_intd%i_Nblocks%i_pot_brush", intsigma, intd, Nblocks); // Don't need a new one of these for every instance.
    brushFile.open(bfilename);
    delete bfilename;

    int counter = 0;
    for(int i=0; i<Nblocks; i++){
        for(int j=0; j<Nblocks; j++){
            centres_x[counter] = (i+1)*intd+(2*i+1)*intsigma + i;
            centres_y[counter] = (j+1)*intd+(2*j+1)*intsigma + j;
            //cout << "centres_x[counter]: " << centres_x[counter] << "; centres_y[counter]: " << centres_y[counter] << endl;
            brushFile << centres_x[counter] << " " << centres_y[counter] << endl;
            counter++;
        }
    }
    brushFile.close();

    //---- Partitioning set up ----//
    vector<int> startpoints = give_startpoints(Nrun, Nsections);
    int len_sections = Nrun/Nsections; // I have already done this. Maybe having the function is not that handy...

    //---- Set up of run ----//
    // Random number distributions
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,Nmat-1);
    std::uniform_int_distribution<int> steps_distr(0,1); // Move: +/-1. Should I be able to keep still too?
    std::uniform_int_distribution<int> dir_distr(0,1);
    std::uniform_real_distribution<double> accept_prob(0,1); // Probability of accepting step with higher energy

    // Variables
    bool free;
    //double thisR2;    // Is this a conflict? Should this be an int?
    int starti, startj, indi, indj, nextindi, nextindj, dir, step, thisR2, nx, ny;
    vector<double> walk_R2(Nrun);                                           // Average R^2
    vector<double> walk_R2_rms(Nrun);                                       // RMS R^2
    vector<vector<double>> walk_R2_store(Nrealizations, vector<double>(Nrun)); // For calculation of the standard deviation
    vector<vector<double>> walk_x(Nrealizations, vector<double>(Nrun));        // Useful in case I want to look at the walks for different starting points (post-processing)
    vector<vector<double>> walk_y(Nrealizations, vector<double>(Nrun));
    vector<vector<double>> walk_R2_sections(Nsections, vector<double>(len_sections));                                                      // Sections
    vector<vector<double>> walk_R2_sections_rms(Nsections, vector<double>(len_sections));                                                  // Sections, rms
    vector<vector<vector<double>>> walk_R2_sections_store(Nrealizations, vector<vector<double>>(Nsections, vector<double>(len_sections))); // For calculation of the standard deviation


    //---- Performing walk ----//
    bool last_locked;
    int consecutive_locked;
    double e_curr, e_next, de, boltzmann, acceptance, illegalmovesuggested;
    nx = 0; ny = 0; acceptance = 0;
    for(int i=0; i<Nrealizations; i++){
        // Resetting some diagnostics quantities:
        last_locked = false;
        consecutive_locked = 0;
        // Start the walk by placing the bead in a random position
        // and check that the position is not occupied
        free = false;
        while(!free){
            indi = distribution(generator); //location();//indi = rand() % Nmat; // Random numbers this way (from <random>) is a bit better than rand() from cstdlib, it seems.
            indj = distribution(generator); //location();//indj = rand() % Nmat;
            //cout << "Initiation. indi: " << indi << "; indj: " << indj << endl;
            e_curr = give_energy(indi, indj, Npots, sigma, power, soften, centres_x, centres_y);
            boltzmann = exp(-beta*e_curr);
            //cout << "walkmat["<<indi<<"]["<<indj<< "]: " << walkmat[indi][indj] << endl;
            //if(walkmat[indi][indj]==0){free=true;} // Only place the bead on a free site
            //cout << "i = " << i << "; Energy: " << e_curr << "; Boltzmann factor:" << boltzmann << endl;
            if(boltzmann>0.5){free=true;} // I've just set some arbitrary threshold so that it won't get too close to a bead.
        }
        starti = indi;
        startj = indj;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl;
        for(int j=0; j<Nrun; j++){
            step = 2*steps_distr(generator)-1;
            //cout << "step: " << step << endl;
            dir = dir_distr(generator);
            //cout << "dir: " << dir << endl;
            if(dir==0){
                //cout << "dir x chosen" << endl;
                nextindi = indi + step;
                nextindj = indj;
                //cout << "nextindi: " << nextindi << "; nextindj: " << nextindj << endl;
                if(nextindi>(Nmat-1) || nextindi<0){
                    nextindi = indi;
                    nextindj = indj;
                    illegalmovesuggested += 1;
                    //cout << "illegalmovesuggested" << endl;
                }
            }
            else{
                //cout << "dir y chosen" << endl;
                nextindi = indi;
                nextindj = indj + step;
                //cout << "nextindi: " << nextindi << "; nextindj: " << nextindj << endl;
                if(nextindj>(Nmat-1) || nextindj<0){
                    nextindi = indi;
                    nextindj = indj;
                    illegalmovesuggested += 1;
                    //cout << "illegalmovesuggested" << endl;
                }
            }
            e_next    = give_energy(nextindi, nextindj, Npots, sigma, power, soften, centres_x, centres_y);
            if(e_next<e_curr){
                // Accept move
                indi = nextindi;
                indj = nextindj;
                e_curr = e_next;
                //cout << "Energy changed, e_curr: " << e_curr << "; e_next: " << e_next << ": e_next lower" << endl;
                acceptance += 1;
                last_locked = false;
                consecutive_locked = 0;
            }
            else{
                de        = e_next-e_curr;
                boltzmann = exp(-beta*de);
                if(accept_prob(generator)<boltzmann){
                    //cout << "Energy changed, e_next: " << e_next << "; e_curr: " << e_curr << ": by prob" << endl;
                    indi = nextindi;
                    indj = nextindj;
                    e_curr = e_next;
                    acceptance += 1;
                    last_locked = false;
                    consecutive_locked = 0;
                    //cout << "Energy changed, e_curr: " << e_curr << endl;
                }
                else{
                    if(last_locked){
                        consecutive_locked++;
                        //cout << "cout: consecutive_locked:" << consecutive_locked << "Energy not changed, e_curr: " << e_curr << "; e_next: " << e_next << "; Current position: [" << indi << "," << indj << "]" << endl;
                        /*
                        if(consecutive_locked>1000){
                            cout << "Energy not changed, e_curr: " << e_curr << "; e_next: " << e_next << "; Current position: [" << indi << "," << indj << "]" << endl;
                        }
                        //*/
                    }
                    /*else{
                        cout << "Move failed, e_next: " << e_next << "; e_curr: " << e_curr << endl;
                    }*/
                    last_locked = true;
                }
            }
            // Average moves
            thisR2 = (indi-starti)*(indi-starti) + (indj-startj)*(indj-startj);
            walk_R2[j] += thisR2;
            walk_R2_store[i][j] = thisR2;
            // Saving coordinates
            walk_x[i][j] = indi;
            walk_y[i][j] = indj;
        }
    }
    acceptance /= (Nrun*Nrealizations);

    for(int i=0; i<Nrun; i++){
        walk_R2[i]/=Nrealizations;
        //cout << "walk_R2[" << i << "]: " << walk_R2[i] << endl;
    } // Averaging complete


    // Finding rms values
    for(int j=0; j<Nrun; j++){
        for(int i=0; i<Nrealizations; i++){
            walk_R2_rms[j] += (walk_R2_store[i][j]-walk_R2[j])*(walk_R2_store[i][j]-walk_R2[j]);
        }
        walk_R2_rms[j] = sqrt(walk_R2_rms[j]/(Nrealizations-1));
    }

    //cout << "nx: " << nx << endl << "ny: " << ny << endl;
    //cout << "Ntotblocks-1: " << Ntotblocks-1 << endl;


    //----- Partitioning -----//

    int startpoint;
    double startx, starty, currx, curry, dx, dy, dotprod;
    for(int i=0; i<Nrealizations; i++){
        for(int j=0; j<Nsections; j++){
            // Extract starting pointsillegalmovesuggestedillegalmovesuggestedillegalmovesuggestedillegalmovesuggestedillegalmovesuggested
            startpoint = startpoints[j];
            startx = walk_x[i][startpoint];
            starty = walk_y[i][startpoint];
            for(int k=1; k<len_sections; k++){ // Start at 1 since first point in array will be 0 anyways
                currx   = walk_x[i][startpoint+k];
                curry   = walk_y[i][startpoint+k];
                dx      = currx-startx;
                dy      = curry-starty;
                dotprod = dx*dx+dy*dy;
                walk_R2_sections[j][k]          += dotprod;
                walk_R2_sections_store[i][j][k]  = dotprod;
            } // End loop over steps (len_sections)
        } // End loop over sections
    } // End loop over realizations

   // Averaging
    for(int i = 0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            walk_R2_sections[i][j] /= Nrealizations;
        } // End loop over steps (len_sections)
    } // End loop over sections

    // rms
    for(int i=0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            for(int k=0; k<Nrealizations; k++){
                walk_R2_sections_rms[i][j] = (walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j])*(walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j]);
            } // End loop over realizations
            walk_R2_sections_rms[i][j] = sqrt(walk_R2_sections_rms[i][j]/(Nrealizations-1));
        } // End loop over steps (len_sections)
    } // End loop over sections


    //----- Writing to file -----//
    //string outfilenamePrefix;


    ofstream outFile;
    char *filename = new char[100000]; // Possibly a long file name
    sprintf(filename, "Nsteps%i_Nreal%i_Nblocks%i_intsigma%i_intd%i_pot_exp%f_sigma%f_factor%f_R2_basic_mc", Nrun, Nrealizations, Nblocks, intsigma, intd, power, sigma, soften);
    outFile.open(filename);
    delete filename;


    for(int i=0; i<Nrun; i++){
        outFile << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl;
        //outFile << std::setprecision(std::numeric_limits<double>::digits10+1) << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl; // Is this long printing really neccessary? // Possibly (probably?) gonna use this later.
    }
    outFile.close();

    ofstream poutFile;
    char *pfilename = new char[100000]; // Possibly a long file name
    sprintf(pfilename, "Nsteps%i_Nreal%i_Nblocks%i_intsigma%i_intd%i_Npart%i_pot_exp%f_sigma%f_factor%f_R2_basic_mc", Nrun, Nrealizations, Nblocks, intsigma, intd, Nsections, power, sigma, soften);
    poutFile.open(pfilename);
    delete pfilename;

    for(int i=0; i<Nsections; i++){
        poutFile << "--- Section " << i <<" ---" << endl;
        for(int j=0; j<len_sections; j++){
            poutFile << j << " " << walk_R2_sections[i][j] << " " << walk_R2_sections_rms[i][j] << endl;
        } // End loop over steps (len_sections)
    } // End loop over sections
    poutFile.close();

    ofstream coordFile;
    char *cfilename = new char[100000]; // Possibly a long file name
    sprintf(cfilename, "Nsteps%i_Nreal%i_Npart%i_pot_exp%f_sigma%f_factor%f_R2_basic_mc_coords", Nrun, Nrealizations, Nsections, power, sigma, soften);
    coordFile.open(cfilename);
    delete cfilename;

    for(int i=0; i<Nrealizations; i++){
        coordFile << "-- Realization " << i << "--" << endl;
        for(int j=0; j<Nrun; j++){
            coordFile << j << " " << walk_x[i][j] << " " << walk_y[i][j] << endl;
        }
    }
    coordFile.close();

    cout << "acceptance rate: " << acceptance << endl;
    cout << "consecutive_locked: " << consecutive_locked << endl;
    cout << "illegalmovesuggested: " << illegalmovesuggested << endl;
    cout << "e_curr: " << e_curr << endl;



}


//-- Periodic BCs --//
void randomwalk_PBC(int sigma, int d, int Nblocks, int Nrun, int Nrealizations, int Nsections, int printevery, bool pbc_codetesting)
{   // Assuming quadratic matrix
    // Sigma: 'radius' of obstactle (i.e. half the length since it is quadratic)
    //---- Matrix set up ----//
    int Nmat, blocksize;
    blocksize = 2*sigma+1;
    Nmat = Nblocks*d + Nblocks*blocksize; // Size of matrix: contribution from spacings + contributions from obstacles (Structure: sp-obst-sp-obst-...-sp-obst)

    //---- Partitioning set up ----//
    vector<int> startpoints = give_startpoints(Nrun, Nsections);
    int len_sections = Nrun/Nsections; // I have already done this. Maybe having the function is not that handy...

    //---- Set up of run ----//
    // Random number distributions
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,Nmat-1);
    std::uniform_int_distribution<int> steps_distr(0,1); // Move: +/-1. Should I be able to keep still too?
    std::uniform_int_distribution<int> dir_distr(0,1);

    // Variables
    bool moveit;
    //double thisR2;    // Is this a conflict? Should this be an int?
    int starti, startj, indi, indj, indi_calc, indj_calc, nextindi, nextindj, dir, step, thisR2, nx, ny, accflagx, accflagy;
    vector<double> walk_R2(Nrun);                                           // Average R^2
    vector<double> walk_R2_rms(Nrun);                                       // RMS R^2
    vector<vector<double>> walk_R2_store(Nrealizations, vector<double>(Nrun)); // For calculation of the standard deviation
    vector<vector<double>> walk_x(Nrealizations, vector<double>(Nrun));        // Useful in case I want to look at the walks for different starting points (post-processing)
    vector<vector<double>> walk_y(Nrealizations, vector<double>(Nrun));
    vector<vector<double>> passing_x(Nrealizations, vector<double>(Nrun));        // Useful in case I want to look at the walks for different starting points (post-processing)
    vector<vector<double>> passing_y(Nrealizations, vector<double>(Nrun));        // To handle boundary crossings
    vector<vector<double>> walk_R2_sections(Nsections, vector<double>(len_sections));                                                      // Sections
    vector<vector<double>> walk_R2_sections_rms(Nsections, vector<double>(len_sections));                                                  // Sections, rms
    vector<vector<vector<double>>> walk_R2_sections_store(Nrealizations, vector<vector<double>>(Nsections, vector<double>(len_sections))); // For calculation of the standard deviation


    //---- Performing walk ----//
    nx = 0; ny = 0;
    for(int i=0; i<Nrealizations; i++){
        // Start the walk by placing the bead in a random position
        // and check that the position is not occupied
        indi = distribution(generator); //location();//indi = rand() % Nmat; // Random numbers this way (from <random>) is a bit better than rand() from cstdlib, it seems.
        indj = distribution(generator); //location();//indj = rand() % Nmat;
        starti = indi;
        startj = indj;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl;
        accflagx = 0;
        accflagy = 0;
        for(int j=0; j<Nrun; j++){
            moveit = false;
            // Just in case:
            passing_x[i][j] = 0;
            passing_y[i][j] = 0;
            while(!moveit){    // Should I let the walker stand still if the move is not accepted?
                step = 2*steps_distr(generator)-1;
                //cout << "step: " << step << endl;
                dir = dir_distr(generator);
                //cout << "dir: " << dir << endl;
                if(dir==0){
                    //cout << "dir x chosen" << endl;
                    nextindi = indi + step;
                    nextindj = indj;
                    moveit=true;
                    nx++;
                    if(nextindi>(Nmat-1)){
                        nextindi  = 0;
                        accflagx += Nmat;
                        passing_x[i][j] = 1;
                    }
                    else if(nextindi<0){
                        nextindi  = Nmat-1;
                        accflagx -=Nmat;
                        passing_x[i][j] = -1;
                    }
                }
                else{
                    //cout << "dir y chosen" << endl;
                    nextindi = indi;
                    nextindj = indj + step;
                    if(nextindj>(Nmat-1)){
                        nextindj  = 0;
                        accflagy += Nmat;
                        passing_y[i][j] = 1;
                    }
                    else if(nextindj<0){
                        nextindj  = Nmat-1;
                        accflagy -= Nmat;
                        passing_y[i][j] = -1;
                    }
                    moveit=true;
                    ny++;
                }

            }
            // Move the walker
            indi = nextindi;
            indj = nextindj;
            // Boundary conditions
            indi_calc = indi + accflagx;
            indj_calc = indj + accflagy;
            // store data somewhere
            // Average moves
            thisR2 = (indi_calc-starti)*(indi_calc-starti) + (indj_calc-startj)*(indj_calc-startj);
            walk_R2[j] += thisR2;
            walk_R2_store[i][j] = thisR2;
            // Saving coordinates
            walk_x[i][j] = indi;
            walk_y[i][j] = indj;
        }
        //cout << "indi: " << indi << ", indj: " << indj << endl;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl << "-------------------------------------" << endl;
    }

    for(int i=0; i<Nrun; i++){
        walk_R2[i]/=Nrealizations;
        //cout << "walk_R2[" << i << "]: " << walk_R2[i] << endl;
    } // Averaging complete


    // Finding rms values
    for(int j=0; j<Nrun; j++){
        for(int i=0; i<Nrealizations; i++){
            walk_R2_rms[j] += (walk_R2_store[i][j]-walk_R2[j])*(walk_R2_store[i][j]-walk_R2[j]);
        }
        walk_R2_rms[j] = sqrt(walk_R2_rms[j]/(Nrealizations-1));
    }

    //cout << "nx: " << nx << endl << "ny: " << ny << endl;
    //cout << "Nmat-1: " << Nmat-1 << endl;


    //----- Partitioning -----//

    int startpoint, thispoint;
    double startx, starty, currx, curry, dx, dy, dotprod;
    for(int i=0; i<Nrealizations; i++){
        for(int j=0; j<Nsections; j++){
            // Extract starting points
            startpoint = startpoints[j];
            startx = walk_x[i][startpoint];
            starty = walk_y[i][startpoint];
            // Reset accflags:
            accflagx = 0;
            accflagy = 0;
            for(int k=1; k<len_sections; k++){ // Start at 1 since first point in array will be 0 anyways
                thispoint = startpoint+k;
                if(passing_x[i][thispoint]!=0){accflagx += passing_x[i][thispoint]*Nmat;} // Add + or - Nmat if a border is crossed
                if(passing_y[i][thispoint]!=0){accflagy += passing_y[i][thispoint]*Nmat;}
                currx   = walk_x[i][startpoint+k]+accflagx;
                curry   = walk_y[i][startpoint+k]+accflagy;
                dx      = currx-startx;
                dy      = curry-starty;
                dotprod = dx*dx+dy*dy;
                walk_R2_sections[j][k]          += dotprod;
                walk_R2_sections_store[i][j][k]  = dotprod;
            } // End loop over steps (len_sections)
        } // End loop over sections
    } // End loop over realizations

   // Averaging
    for(int i = 0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            walk_R2_sections[i][j] /= Nrealizations;
        } // End loop over steps (len_sections)
    } // End loop over sections

    // rms
    for(int i=0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            for(int k=0; k<Nrealizations; k++){
                walk_R2_sections_rms[i][j] = (walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j])*(walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j]);
            } // End loop over realizations
            walk_R2_sections_rms[i][j] = sqrt(walk_R2_sections_rms[i][j]/(Nrealizations-1));
        } // End loop over steps (len_sections)
    } // End loop over sections


    //----- Writing to file -----//
    //string outfilenamePrefix;


    ofstream outFile;
    char *filename = new char[100000]; // Possibly a long file name
    sprintf(filename, "PBC/sigma%i_d%i/Nblocks%i/randomwalk_R2_Nsteps%i_Nreal%i_printevery%i", sigma, d, Nblocks, Nrun, Nrealizations,printevery);//Nrun, Nrealizations, sigma, d, Nblocks);
    //sprintf(filename, "Nsteps%i_Nreal%i_sigma%i_d%i_Nblocks%i_randomwalk_PBC_R2_basic_mc", Nrun, Nrealizations, sigma, d, Nblocks);
    outFile.open(filename);
    delete filename;


    for(int i=0; i<Nrun; i++){
        if(i%printevery==0){
            outFile << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl;
        }
        //outFile << std::setprecision(std::numeric_limits<double>::digits10+1) << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl; // Is this long printing really neccessary? // Possibly (probably?) gonna use this later.
    }
    outFile.close();

    ofstream poutFile;
    char *pfilename = new char[100000]; // Possibly a long file name
    sprintf(pfilename, "PBC/sigma%i_d%i/Nblocks%i/randomwalk_R2_Nsteps%i_Nreal%i_Npart%i_printevery%i", sigma, d, Nblocks, Nrun, Nrealizations, Nsections,printevery);//Nrun, Nrealizations, sigma, d, Nblocks);
    //sprintf(pfilename, "Nsteps%i_Nreal%i_sigma%i_d%i_Nblocks%i_Npart%i_randomwalk_PBC_R2_basic_mc", Nrun, Nrealizations, sigma, d, Nblocks, Nsections);
    poutFile.open(pfilename);
    delete pfilename;

    for(int i=0; i<Nsections; i++){
        poutFile << "--- Section " << i <<" ---" << endl;
        for(int j=0; j<len_sections; j++){
            if(j%printevery==0){
                poutFile << j << " " << walk_R2_sections[i][j] << " " << walk_R2_sections_rms[i][j] << endl;
            }
        } // End loop over steps (len_sections)
    } // End loop over sections
    poutFile.close();

    if(pbc_codetesting){
        ofstream coordFile;
        char *cfilename = new char[100000]; // Possibly a long file name
        sprintf(cfilename, "PBC/sigma%i_d%i/Nblocks%i/randomwalk_Nsteps%i_Nreal%i_coords", sigma, d, Nblocks, Nrun, Nrealizations);
        //sprintf(cfilename, "Nsteps%i_Nreal%i_Npart%i_pot_exp%f_sigma%f_factor%f_PBC_R2_basic_mc_coords", Nrun, Nrealizations, Nsections, power, sigma, soften);
        coordFile.open(cfilename);
        delete cfilename;

        for(int i=0; i<Nrealizations; i++){
            coordFile << "-- Realization " << i << "--" << endl;
            for(int j=0; j<Nrun; j++){
                coordFile << j << " " << walk_x[i][j] << " " << walk_y[i][j] << " ; " << passing_x[i][j] << " " << passing_y[i][j] << endl;
            }
        }
        coordFile.close();
    }
}



void matrixwalk_hard_easyinput_PBC(int sigma, int d, int Nblocks, int Nrun, int Nrealizations, int Nsections, int printevery, bool printmatrix, bool pbc_codetesting)
{   // Assuming quadratic matrix
    // Sigma: 'radius' of obstactle (i.e. half the length since it is quadratic)
    //---- Matrix set up ----//
    int Nmat, blocksize, spacingsi, spacingsj;
    blocksize = 2*sigma+1;
    Nmat = Nblocks*d + Nblocks*blocksize; // Size of matrix: contribution from spacings + contributions from obstacles (Structure: sp-obst-sp-obst-...-sp-obst)
    vector<vector<int>> walkmat(Nmat, vector<int>(Nmat));

    // Setting matrix blocks
    // Must be a handier way to do this?:
    for(int i=0; i<Nblocks; i++){
        for(int j=0; j<Nblocks; j++){
            for(int k=0; k<blocksize; k++){
                for(int l=0; l<blocksize; l++){
                    spacingsi = (i+0.5)*d; // Have to make sure that index is an int
                    spacingsj = (j+0.5)*d; // Have to make sure that index is an int
                    walkmat[spacingsi+i*blocksize+k][spacingsj+j*blocksize+l] = 1;
                }
            }
        }
    }

    //---- Partitioning set up ----//
    vector<int> startpoints = give_startpoints(Nrun, Nsections);
    int len_sections = Nrun/Nsections; // I have already done this. Maybe having the function is not that handy...

    //---- Set up of run ----//
    // Random number distributions
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,Nmat-1);
    std::uniform_int_distribution<int> steps_distr(0,1); // Move: +/-1. Should I be able to keep still too?
    std::uniform_int_distribution<int> dir_distr(0,1);

    // Variables
    bool free, moveit;
    //double thisR2;    // Is this a conflict? Should this be an int?
    int starti, startj, indi, indj, indi_calc, indj_calc, nextindi, nextindj, dir, step, thisR2, nx, ny, accflagx, accflagy;
    vector<double> walk_R2(Nrun);                                           // Average R^2
    vector<double> walk_R2_rms(Nrun);                                       // RMS R^2
    vector<vector<double>> walk_R2_store(Nrealizations, vector<double>(Nrun)); // For calculation of the standard deviation
    vector<vector<double>> walk_x(Nrealizations, vector<double>(Nrun));        // Useful in case I want to look at the walks for different starting points (post-processing)
    vector<vector<double>> walk_y(Nrealizations, vector<double>(Nrun));
    vector<vector<double>> passing_x(Nrealizations, vector<double>(Nrun));        // Useful in case I want to look at the walks for different starting points (post-processing)
    vector<vector<double>> passing_y(Nrealizations, vector<double>(Nrun));        // To handle boundary crossings
    vector<vector<double>> walk_R2_sections(Nsections, vector<double>(len_sections));                                                      // Sections
    vector<vector<double>> walk_R2_sections_rms(Nsections, vector<double>(len_sections));                                                  // Sections, rms
    vector<vector<vector<double>>> walk_R2_sections_store(Nrealizations, vector<vector<double>>(Nsections, vector<double>(len_sections))); // For calculation of the standard deviation


    //---- Performing walk ----//
    nx = 0; ny = 0;
    for(int i=0; i<Nrealizations; i++){
        // Start the walk by placing the bead in a random position
        // and check that the position is not occupied
        free = false;
        while(!free){
            indi = distribution(generator); //location();//indi = rand() % Nmat; // Random numbers this way (from <random>) is a bit better than rand() from cstdlib, it seems.
            indj = distribution(generator); //location();//indj = rand() % Nmat;
            //cout << "walkmat["<<indi<<"]["<<indj<< "]: " << walkmat[indi][indj] << endl;
            if(walkmat[indi][indj]==0){free=true;} // Only place the bead on a free site
        }
        starti = indi;
        startj = indj;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl;
        accflagx = 0;
        accflagy = 0;
        for(int j=0; j<Nrun; j++){
            // Just in case:
            passing_x[i][j] = 0;
            passing_y[i][j] = 0;
            moveit = false;
            while(!moveit){    // Should I let the walker stand still if the move is not accepted?
                step = 2*steps_distr(generator)-1;
                //cout << "step: " << step << endl;
                dir = dir_distr(generator);
                //cout << "dir: " << dir << endl;
                if(dir==0){
                    //cout << "dir x chosen" << endl;
                    nextindi = indi + step;
                    nextindj = indj;
                    if(nextindi>(Nmat-1)){
                        nextindi = 0; // No test neccessary. The sites at the edges are filled.
                        accflagx += Nmat;
                        passing_x[i][j] = 1;
                        moveit=true;
                        nx++;
                    }
                    else if(nextindi<0){
                        nextindi = Nmat-1;
                        accflagx -= Nmat;
                        passing_x[i][j] = -1;
                        moveit=true;
                        nx++;
                    }
                    else{  //if(nextindi<Nmat && nextindi>-1){ // I shouldn't need this test because of the other ones.
                        if(walkmat[nextindi][nextindj]==0){
                            moveit=true;
                            nx++;
                        }
                    }
                }
                else{
                    //cout << "dir y chosen" << endl;
                    nextindi = indi;
                    nextindj = indj + step;
                    if(nextindj>(Nmat-1)){
                        nextindj = 0; // No test neccessary. The sites at the edges are filled.
                        accflagy += Nmat;
                        passing_y[i][j] = 1;
                        moveit=true;
                        ny++;
                    }
                    else if(nextindj<0){
                        nextindj = Nmat-1;
                        accflagy -= Nmat;
                        passing_y[i][j] = -1;
                        moveit=true;
                        ny++;
                    }
                    else{  //if(nextindj<Nmat && nextindj>-1){ // I shouldn't need this test because of the other ones.
                        if(walkmat[nextindi][nextindj]==0){
                            moveit=true;
                            ny++;
                        }
                    }
                }
            }
            // Move the walker
            indi = nextindi;
            indj = nextindj;
            // Boundary conditions
            indi_calc = indi + accflagx;
            indj_calc = indj + accflagy;
            // store data somewhere
            // Average moves
            thisR2 = (indi_calc-starti)*(indi_calc-starti) + (indj_calc-startj)*(indj_calc-startj);
            walk_R2[j] += thisR2;
            walk_R2_store[i][j] = thisR2;
            // Saving coordinates
            walk_x[i][j] = indi;
            walk_y[i][j] = indj;
        }
        //cout << "indi: " << indi << ", indj: " << indj << endl;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl << "-------------------------------------" << endl;
    }

    for(int i=0; i<Nrun; i++){
        walk_R2[i]/=Nrealizations;
        //cout << "walk_R2[" << i << "]: " << walk_R2[i] << endl;
    } // Averaging complete


    // Finding rms values
    for(int j=0; j<Nrun; j++){
        for(int i=0; i<Nrealizations; i++){
            walk_R2_rms[j] += (walk_R2_store[i][j]-walk_R2[j])*(walk_R2_store[i][j]-walk_R2[j]);
        }
        walk_R2_rms[j] = sqrt(walk_R2_rms[j]/(Nrealizations-1));
    }

    //cout << "nx: " << nx << endl << "ny: " << ny << endl;
    //cout << "Nmat-1: " << Nmat-1 << endl;


    //----- Partitioning -----//

    int startpoint, thispoint;
    double startx, starty, currx, curry, dx, dy, dotprod;
    for(int i=0; i<Nrealizations; i++){
        for(int j=0; j<Nsections; j++){
            // Extract starting points
            startpoint = startpoints[j];
            startx = walk_x[i][startpoint];
            starty = walk_y[i][startpoint];
            // Reset accflags:
            accflagx = 0;
            accflagy = 0;
            for(int k=1; k<len_sections; k++){ // Start at 1 since first point in array will be 0 anyways
                thispoint = startpoint+k;
                if(passing_x[i][thispoint]!=0){accflagx += passing_x[i][thispoint]*Nmat;} // Add + or - Nmat if a border is crossed
                if(passing_y[i][thispoint]!=0){accflagy += passing_y[i][thispoint]*Nmat;}
                currx   = walk_x[i][startpoint+k]+accflagx;
                curry   = walk_y[i][startpoint+k]+accflagy;
                dx      = currx-startx;
                dy      = curry-starty;
                dotprod = dx*dx+dy*dy;
                walk_R2_sections[j][k]          += dotprod;
                walk_R2_sections_store[i][j][k]  = dotprod;
            } // End loop over steps (len_sections)
        } // End loop over sections
    } // End loop over realizations

   // Averaging
    for(int i = 0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            walk_R2_sections[i][j] /= Nrealizations;
        } // End loop over steps (len_sections)
    } // End loop over sections

    // rms
    for(int i=0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            for(int k=0; k<Nrealizations; k++){
                walk_R2_sections_rms[i][j] = (walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j])*(walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j]);
            } // End loop over realizations
            walk_R2_sections_rms[i][j] = sqrt(walk_R2_sections_rms[i][j]/(Nrealizations-1));
        } // End loop over steps (len_sections)
    } // End loop over sections


    //----- Writing to file -----//
    //string outfilenamePrefix;


    ofstream outFile;
    char *filename = new char[100000]; // Possibly a long file name
    sprintf(filename, "PBC/sigma%i_d%i/Nblocks%i/hardpotwalk_R2_Nsteps%i_Nreal%i_printevery%i", sigma, d, Nblocks, Nrun, Nrealizations, printevery);//Nrun, Nrealizations, sigma, d, Nblocks);
    //sprintf(filename, "Nsteps%i_Nreal%i_sigma%i_d%i_Nblocks%i_hardpot_PBC_R2_basic_mc", Nrun, Nrealizations, sigma, d, Nblocks);
    outFile.open(filename);
    delete filename;


    for(int i=0; i<Nrun; i++){
        if(i%printevery==0){
            outFile << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl;
        }
        //outFile << std::setprecision(std::numeric_limits<double>::digits10+1) << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl; // Is this long printing really neccessary? // Possibly (probably?) gonna use this later.
    }
    outFile.close();

    ofstream poutFile;
    char *pfilename = new char[100000]; // Possibly a long file name
    sprintf(pfilename, "PBC/sigma%i_d%i/Nblocks%i/hardpotwalk_R2_Nsteps%i_Nreal%i_Npart%i_printevery%i", sigma, d, Nblocks, Nrun, Nrealizations, Nsections,printevery);//Nrun, Nrealizations, sigma, d, Nblocks);
    //sprintf(pfilename, "Nsteps%i_Nreal%i_sigma%i_d%i_Nblocks%i_Npart%i_hardpot_PBC_R2_basic_mc", Nrun, Nrealizations, sigma, d, Nblocks, Nsections);
    poutFile.open(pfilename);
    delete pfilename;

    for(int i=0; i<Nsections; i++){
        poutFile << "--- Section " << i <<" ---" << endl;
        for(int j=0; j<len_sections; j++){
            if(j%printevery==0){
                poutFile << j << " " << walk_R2_sections[i][j] << " " << walk_R2_sections_rms[i][j] << endl;
            }
        } // End loop over steps (len_sections)
    } // End loop over sections
    poutFile.close();

    if(printmatrix){
        // Print matrix for testing:
        ofstream moutFile;
        char *mfilename = new char[100000]; // Possibly a long file name
        sprintf(mfilename, "PBC/sigma%i_d%i/Nblocks%i/Matrices_etc/matrix_walk", sigma, d, Nblocks);
        moutFile.open(mfilename);
        delete mfilename;

        for(int i=0; i<Nmat; i++){
            for(int j=0; j<Nmat; j++){
                moutFile << walkmat[i][j] << " ";
            } // End loop over steps (len_sections)
            moutFile << endl;
        } // End loop over sections
        moutFile.close();
    } // End if(printmatrix)

    if(pbc_codetesting){
        ofstream coordFile;
        char *cfilename = new char[100000]; // Possibly a long file name
        sprintf(cfilename, "PBC/sigma%i_d%i/Nblocks%i/matrixwalk_Nsteps%i_Nreal%i_coords", sigma, d, Nblocks, Nrun, Nrealizations);
        coordFile.open(cfilename);
        delete cfilename;

        for(int i=0; i<Nrealizations; i++){
            coordFile << "-- Realization " << i << "--" << endl;
            for(int j=0; j<Nrun; j++){
                coordFile << j << " " << walk_x[i][j] << " " << walk_y[i][j] << " ; " << passing_x[i][j] << " " << passing_y[i][j] << endl;
            }
        }
        coordFile.close();
    }

}

// MC
void matrixmc_hard_easyinput_PBC(int sigma, int d, int Nblocks, int Nrun, int Nrealizations, int Nsections, int printevery, bool printmatrix, bool pbc_codetesting)
{   // Assuming quadratic matrix
    // Sigma: 'radius' of obstactle (i.e. half the length since it is quadratic)
    //---- Matrix set up ----//
    int Nmat, blocksize, spacingsi, spacingsj;
    blocksize = 2*sigma+1;
    Nmat = Nblocks*d + Nblocks*blocksize; // Size of matrix: contribution from spacings + contributions from obstacles (Structure: sp-obst-sp-obst-...-sp-obst)
    vector<vector<int>> walkmat(Nmat, vector<int>(Nmat));

    // Setting matrix blocks
    // Must be a handier way to do this?:
    for(int i=0; i<Nblocks; i++){
        for(int j=0; j<Nblocks; j++){
            for(int k=0; k<blocksize; k++){
                for(int l=0; l<blocksize; l++){
                    spacingsi = (i+0.5)*d; // Have to make sure that index is an int
                    spacingsj = (j+0.5)*d; // Have to make sure that index is an int
                    walkmat[spacingsi+i*blocksize+k][spacingsj+j*blocksize+l] = 1;
                }
            }
        }
    }

    //---- Partitioning set up ----//
    vector<int> startpoints = give_startpoints(Nrun, Nsections);
    int len_sections = Nrun/Nsections; // I have already done this. Maybe having the function is not that handy...

    //---- Set up of run ----//
    // Random number distributions
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,Nmat-1);
    std::uniform_int_distribution<int> steps_distr(0,1); // Move: +/-1. Should I be able to keep still too?
    std::uniform_int_distribution<int> dir_distr(0,1);

    // Variables
    bool free;
    //double thisR2;    // Is this a conflict? Should this be an int?
    int starti, startj, indi, indj, indi_calc, indj_calc, nextindi, nextindj, dir, step, thisR2, nx, ny, accflagx, accflagy;
    vector<double> walk_R2(Nrun);                                           // Average R^2
    vector<double> walk_R2_rms(Nrun);                                       // RMS R^2
    vector<vector<double>> walk_R2_store(Nrealizations, vector<double>(Nrun)); // For calculation of the standard deviation
    vector<vector<double>> walk_x(Nrealizations, vector<double>(Nrun));        // Useful in case I want to look at the walks for different starting points (post-processing)
    vector<vector<double>> walk_y(Nrealizations, vector<double>(Nrun));
    vector<vector<double>> passing_x(Nrealizations, vector<double>(Nrun));        // Useful in case I want to look at the walks for different starting points (post-processing)
    vector<vector<double>> passing_y(Nrealizations, vector<double>(Nrun));        // To handle boundary crossings
    vector<vector<double>> walk_R2_sections(Nsections, vector<double>(len_sections));                                                      // Sections
    vector<vector<double>> walk_R2_sections_rms(Nsections, vector<double>(len_sections));                                                  // Sections, rms
    vector<vector<vector<double>>> walk_R2_sections_store(Nrealizations, vector<vector<double>>(Nsections, vector<double>(len_sections))); // For calculation of the standard deviation


    //---- Performing walk ----//
    nx = 0; ny = 0;
    for(int i=0; i<Nrealizations; i++){
        // Start the walk by placing the bead in a random position
        // and check that the position is not occupied
        free = false;
        while(!free){
            indi = distribution(generator); //location();//indi = rand() % Nmat; // Random numbers this way (from <random>) is a bit better than rand() from cstdlib, it seems.
            indj = distribution(generator); //location();//indj = rand() % Nmat;
            //cout << "walkmat["<<indi<<"]["<<indj<< "]: " << walkmat[indi][indj] << endl;
            if(walkmat[indi][indj]==0){free=true;} // Only place the bead on a free site
        }
        starti = indi;
        startj = indj;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl;
        accflagx = 0;
        accflagy = 0;
        for(int j=0; j<Nrun; j++){
            // Just in case:
            passing_x[i][j] = 0;
            passing_y[i][j] = 0;
            step = 2*steps_distr(generator)-1;
            //cout << "step: " << step << endl;
            dir = dir_distr(generator);
            //cout << "dir: " << dir << endl;
            if(dir==0){
                //cout << "dir x chosen" << endl;
                nextindi = indi + step;
                nextindj = indj;
                if(nextindi>(Nmat-1)){
                    nextindi = 0; // No test neccessary. The sites at the edges are filled.
                    accflagx += Nmat;
                    passing_x[i][j] = 1;
                    nx++;
                }
                else if(nextindi<0){
                    nextindi = Nmat-1;
                    accflagx -= Nmat;
                    passing_x[i][j] = -1;
                    nx++;
                }
                else{  //if(nextindi<Nmat && nextindi>-1){ // I shouldn't need this test because of the other ones.
                    if(walkmat[nextindi][nextindj]==0){
                        nx++;
                    }
                    else{
                        nextindi = indi;
                        nextindj = indj;
                    }
                }
            }
            else{
                //cout << "dir y chosen" << endl;
                nextindi = indi;
                nextindj = indj + step;
                if(nextindj>(Nmat-1)){
                    nextindj = 0; // No test neccessary. The sites at the edges are filled.
                    accflagy += Nmat;
                    passing_y[i][j] = 1;
                    ny++;
                }
                else if(nextindj<0){
                    nextindj = Nmat-1;
                    accflagy -= Nmat;
                    passing_y[i][j] = -1;
                    ny++;
                }
                else{  //if(nextindi<Nmat && nextindi>-1){ // I shouldn't need this test because of the other ones.
                    if(walkmat[nextindi][nextindj]==0){
                        ny++;
                    }
                    else{
                        nextindi = indi;
                        nextindj = indj;
                    }
                }
            }
            // Move the walker
            indi = nextindi;
            indj = nextindj;
            // Boundary conditions
            indi_calc = indi + accflagx;
            indj_calc = indj + accflagy;
            // store data somewhere
            // Average moves
            thisR2 = (indi_calc-starti)*(indi_calc-starti) + (indj_calc-startj)*(indj_calc-startj);
            walk_R2[j] += thisR2;
            walk_R2_store[i][j] = thisR2;
            // Saving coordinates
            walk_x[i][j] = indi;
            walk_y[i][j] = indj;
        }
        //cout << "indi: " << indi << ", indj: " << indj << endl;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl << "-------------------------------------" << endl;
    }

    for(int i=0; i<Nrun; i++){
        walk_R2[i]/=Nrealizations;
        //cout << "walk_R2[" << i << "]: " << walk_R2[i] << endl;
    } // Averaging complete


    // Finding rms values
    for(int j=0; j<Nrun; j++){
        for(int i=0; i<Nrealizations; i++){
            walk_R2_rms[j] += (walk_R2_store[i][j]-walk_R2[j])*(walk_R2_store[i][j]-walk_R2[j]);
        }
        walk_R2_rms[j] = sqrt(walk_R2_rms[j]/(Nrealizations-1));
    }

    //cout << "nx: " << nx << endl << "ny: " << ny << endl;
    //cout << "Nmat-1: " << Nmat-1 << endl;


    //----- Partitioning -----//

    int startpoint, thispoint;
    double startx, starty, currx, curry, dx, dy, dotprod;
    for(int i=0; i<Nrealizations; i++){
        for(int j=0; j<Nsections; j++){
            // Extract starting points
            startpoint = startpoints[j];
            startx = walk_x[i][startpoint];
            starty = walk_y[i][startpoint];
            // Reset accflags:
            accflagx = 0;
            accflagy = 0;
            for(int k=1; k<len_sections; k++){ // Start at 1 since first point in array will be 0 anyways
                thispoint = startpoint+k;
                if(passing_x[i][thispoint]!=0){accflagx += passing_x[i][thispoint]*Nmat;} // Add + or - Nmat if a border is crossed
                if(passing_y[i][thispoint]!=0){accflagy += passing_y[i][thispoint]*Nmat;}
                currx   = walk_x[i][startpoint+k]+accflagx;
                curry   = walk_y[i][startpoint+k]+accflagy;
                dx      = currx-startx;
                dy      = curry-starty;
                dotprod = dx*dx+dy*dy;
                walk_R2_sections[j][k]          += dotprod;
                walk_R2_sections_store[i][j][k]  = dotprod;
            } // End loop over steps (len_sections)
        } // End loop over sections
    } // End loop over realizations

   // Averaging
    for(int i = 0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            walk_R2_sections[i][j] /= Nrealizations;
        } // End loop over steps (len_sections)
    } // End loop over sections

    // rms
    for(int i=0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            for(int k=0; k<Nrealizations; k++){
                walk_R2_sections_rms[i][j] = (walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j])*(walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j]);
            } // End loop over realizations
            walk_R2_sections_rms[i][j] = sqrt(walk_R2_sections_rms[i][j]/(Nrealizations-1));
        } // End loop over steps (len_sections)
    } // End loop over sections


    //----- Writing to file -----//
    //string outfilenamePrefix;


    ofstream outFile;
    char *filename = new char[100000]; // Possibly a long file name
    sprintf(filename, "PBC/sigma%i_d%i/Nblocks%i/hardpotmc_R2_Nsteps%i_Nreal%i_printevery%i", sigma, d, Nblocks, Nrun, Nrealizations, printevery);//Nrun, Nrealizations, sigma, d, Nblocks);
    //sprintf(filename, "Nsteps%i_Nreal%i_sigma%i_d%i_Nblocks%i_hardpot_PBC_R2_basic", Nrun, Nrealizations, sigma, d, Nblocks);
    outFile.open(filename);
    delete filename;


    for(int i=0; i<Nrun; i++){
        if(i%printevery==0){
           outFile << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl;
        }
        //outFile << std::setprecision(std::numeric_limits<double>::digits10+1) << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl; // Is this long printing really neccessary? // Possibly (probably?) gonna use this later.
    }
    outFile.close();

    ofstream poutFile;
    char *pfilename = new char[100000]; // Possibly a long file name
    sprintf(pfilename, "PBC/sigma%i_d%i/Nblocks%i/hardpotmc_R2_Nsteps%i_Nreal%i_Npart%i_printevery%i", sigma, d, Nblocks, Nrun, Nrealizations, Nsections, printevery);//Nrun, Nrealizations, sigma, d, Nblocks);
    //sprintf(pfilename, "Nsteps%i_Nreal%i_sigma%i_d%i_Nblocks%i_Npart%i_hardpot_PBC_R2_basic", Nrun, Nrealizations, sigma, d, Nblocks, Nsections);
    poutFile.open(pfilename);
    delete pfilename;

    for(int i=0; i<Nsections; i++){
        poutFile << "--- Section " << i <<" ---" << endl;
        for(int j=0; j<len_sections; j++){
            if(j%printevery==0){
                poutFile << j << " " << walk_R2_sections[i][j] << " " << walk_R2_sections_rms[i][j] << endl;
            }
        } // End loop over steps (len_sections)
    } // End loop over sections
    poutFile.close();

    if(printmatrix){
        // Print matrix for testing:
        ofstream moutFile;
        char *mfilename = new char[100000]; // Possibly a long file name
        sprintf(mfilename, "PBC/sigma%i_d%i//Nblocks%i/Matrices_etc/matrix_mc", sigma, d, Nblocks);
        moutFile.open(mfilename);
        delete mfilename;

        for(int i=0; i<Nmat; i++){
            for(int j=0; j<Nmat; j++){
                moutFile << walkmat[i][j] << " ";
            } // End loop over steps (len_sections)
            moutFile << endl;
        } // End loop over sections
        moutFile.close();
    } // End if(printmatrix)

    if(pbc_codetesting){
        ofstream coordFile;
        char *cfilename = new char[100000]; // Possibly a long file name
        sprintf(cfilename, "PBC/sigma%i_d%i/Nblocks%i/matrixmc_Nsteps%i_Nreal%i_coords", sigma, d, Nblocks, Nrun, Nrealizations);
        coordFile.open(cfilename);
        delete cfilename;

        for(int i=0; i<Nrealizations; i++){
            coordFile << "-- Realization " << i << "--" << endl;
            for(int j=0; j<Nrun; j++){
                coordFile << j << " " << walk_x[i][j] << " " << walk_y[i][j] << " ; " << passing_x[i][j] << " " << passing_y[i][j] << endl;
            }
        }
        coordFile.close();
    }

}


void matrixmc_potential_easyinput_PBC(int intsigma, int intd, int Nblocks, int Nrun, int Nrealizations, int Nsections, int printevery, double beta, double sigma, double power, double soften, bool pbc_codetesting)
{
    int Nmat, blocksize, spacingsi, spacingsj;
    blocksize = 2*intsigma+1;
    Nmat = Nblocks*intd + Nblocks*blocksize ; // Size of matrix: contribution from spacings + contributions from obstacles (Structure: sp-obst-sp-obst-...-sp-obst)
    //vector<vector<int>> walkmat(Nmat, vector<int>(Nmat)); // I don't really need walkmat, at least not yet.

    //---- Potential set up ----//
    // Number of 'obstacles' given by Nblocks x Nblocks
    // Centre them evenly in a grid-like manner
    int Npots = Nblocks*Nblocks;
    vector<double> centres_x = vector<double>(Npots);
    vector<double> centres_y = vector<double>(Npots);


    ofstream brushFile;
    char *bfilename = new char[100000]; // Possibly a long file name
    sprintf(bfilename, "PBC/sigma%i_d%i/Nblocks%i/Matrices_etc/pot_brush", intsigma, intd, Nblocks);
    //sprintf(bfilename, "intsigma%i_intd%i_Nblocks%i_pot_brush", intsigma, intd, Nblocks); // Don't need a new one of these for every instance.
    brushFile.open(bfilename);
    delete bfilename;

    int counter = 0;
    for(int i=0; i<Nblocks; i++){
        for(int j=0; j<Nblocks; j++){
            spacingsi = (i+0.5)*intd;
            spacingsj = (j+0.5)*intd;
            centres_x[counter] = spacingsi+(2*i+1)*intsigma + i;
            centres_y[counter] = spacingsj+(2*j+1)*intsigma + j;
            cout << "centres_x[counter]: " << centres_x[counter] << "; centres_y[counter]: " << centres_y[counter] << endl;
            brushFile << centres_x[counter] << " " << centres_y[counter] << endl;
            counter++;
        }
    }
    brushFile.close();

    //---- Partitioning set up ----//
    vector<int> startpoints = give_startpoints(Nrun, Nsections);
    int len_sections = Nrun/Nsections; // I have already done this. Maybe having the function is not that handy...

    //---- Set up of run ----//
    // Random number distributions
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,Nmat-1);
    std::uniform_int_distribution<int> steps_distr(0,1); // Move: +/-1. Should I be able to keep still too?
    std::uniform_int_distribution<int> dir_distr(0,1);
    std::uniform_real_distribution<double> accept_prob(0,1); // Probability of accepting step with higher energy

    // Variables
    bool free;
    //double thisR2;    // Is this a conflict? Should this be an int?
    int starti, startj, indi, indj, indi_calc, indj_calc, nextindi, nextindj, dir, step, thisR2, nx, ny, accflagx, accflagy;
    vector<double> walk_R2(Nrun);                                           // Average R^2
    vector<double> walk_R2_rms(Nrun);                                       // RMS R^2
    vector<vector<double>> walk_R2_store(Nrealizations, vector<double>(Nrun)); // For calculation of the standard deviation
    vector<vector<double>> walk_x(Nrealizations, vector<double>(Nrun));        // Useful in case I want to look at the walks for different starting points (post-processing)
    vector<vector<double>> walk_y(Nrealizations, vector<double>(Nrun));
    vector<vector<double>> passing_x(Nrealizations, vector<double>(Nrun));        // Useful in case I want to look at the walks for different starting points (post-processing)
    vector<vector<double>> passing_y(Nrealizations, vector<double>(Nrun));        // To handle boundary crossings
    vector<vector<double>> walk_R2_sections(Nsections, vector<double>(len_sections));                                                      // Sections
    vector<vector<double>> walk_R2_sections_rms(Nsections, vector<double>(len_sections));                                                  // Sections, rms
    vector<vector<vector<double>>> walk_R2_sections_store(Nrealizations, vector<vector<double>>(Nsections, vector<double>(len_sections))); // For calculation of the standard deviation


    //---- Performing walk ----//
    bool last_locked;
    int consecutive_locked;
    double e_curr, e_next, de, boltzmann, acceptance;
    nx = 0; ny = 0; acceptance = 0;
    for(int i=0; i<Nrealizations; i++){
        // Resetting some diagnostics quantities:
        last_locked = false;
        consecutive_locked = 0;
        // Start the walk by placing the bead in a random position
        // and check that the position is not occupied
        free = false;
        while(!free){
            indi = distribution(generator); //location();//indi = rand() % Nmat; // Random numbers this way (from <random>) is a bit better than rand() from cstdlib, it seems.
            indj = distribution(generator); //location();//indj = rand() % Nmat;
            //cout << "Initiation. indi: " << indi << "; indj: " << indj << endl;
            e_curr = give_energy(indi, indj, Npots, sigma, power, soften, centres_x, centres_y);
            boltzmann = exp(-beta*e_curr);
            //cout << "walkmat["<<indi<<"]["<<indj<< "]: " << walkmat[indi][indj] << endl;
            //if(walkmat[indi][indj]==0){free=true;} // Only place the bead on a free site
            //cout << "i = " << i << "; Energy: " << e_curr << "; Boltzmann factor:" << boltzmann << endl;
            if(boltzmann>0.5){free=true;} // I've just set some arbitrary threshold so that it won't get too close to a bead.
        }
        starti = indi;
        startj = indj;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl;
        accflagx = 0;
        accflagy = 0;
        for(int j=0; j<Nrun; j++){
            // Just in case:
            passing_x[i][j] = 0;
            passing_y[i][j] = 0;
            step = 2*steps_distr(generator)-1;
            //cout << "step: " << step << endl;
            dir = dir_distr(generator);
            //cout << "dir: " << dir << endl;
            if(dir==0){
                //cout << "dir x chosen" << endl;
                nextindi = indi + step;
                nextindj = indj;
                //cout << "nextindi: " << nextindi << "; nextindj: " << nextindj << endl;
                if(nextindi>Nmat-1){
                    nextindi = 0;
                    accflagx += Nmat;
                    passing_x[i][j] = 1;
                }
                else if(nextindi<0){
                    nextindi = Nmat-1;
                    accflagx -= Nmat;
                    passing_x[i][j] = -1;
                }
            }
            else{
                //cout << "dir y chosen" << endl;
                nextindi = indi;
                nextindj = indj + step;
                //cout << "nextindi: " << nextindi << "; nextindj: " << nextindj << endl;
                if(nextindj>Nmat-1){
                    nextindj = 0;
                    accflagy += Nmat;
                    passing_y[i][j] = 1;
                }
                else if(nextindj<0){
                    nextindj = Nmat-1;
                    accflagy -= Nmat;
                    passing_y[i][j] = -1;
                }
            }
            e_next    = give_energy(nextindi, nextindj, Npots, sigma, power, soften, centres_x, centres_y);
            if(e_next<e_curr){
                // Accept move
                indi = nextindi;
                indj = nextindj;
                e_curr = e_next;
                //cout << "Energy changed, e_curr: " << e_curr << "; e_next: " << e_next << ": e_next lower" << endl;
                acceptance += 1;
                last_locked = false;
                consecutive_locked = 0;
            }
            else{
                de        = e_next-e_curr;
                boltzmann = exp(-beta*de);
                if(accept_prob(generator)<boltzmann){
                    //cout << "Energy changed, e_next: " << e_next << "; e_curr: " << e_curr << ": by prob" << endl;
                    indi = nextindi;
                    indj = nextindj;
                    e_curr = e_next;
                    acceptance += 1;
                    last_locked = false;
                    consecutive_locked = 0;
                    //cout << "Energy changed, e_curr: " << e_curr << endl;
                }
                else{
                    if(last_locked){
                        consecutive_locked++;
                        //cout << "cout: consecutive_locked:" << consecutive_locked << "Energy not changed, e_curr: " << e_curr << "; e_next: " << e_next << "; Current position: [" << indi << "," << indj << "]" << endl;
                        /*
                        if(consecutive_locked>1000){
                            cout << "Energy not changed, e_curr: " << e_curr << "; e_next: " << e_next << "; Current position: [" << indi << "," << indj << "]" << endl;
                        }
                        //*/
                    }
                    /*else{
                        cout << "Move failed, e_next: " << e_next << "; e_curr: " << e_curr << endl;
                    }*/
                    last_locked = true;
                }
            }
            // Boundary conditions
            indi_calc = indi + accflagx;
            indj_calc = indj + accflagy;
            // store data somewhere
            // Average moves
            thisR2 = (indi_calc-starti)*(indi_calc-starti) + (indj_calc-startj)*(indj_calc-startj);
            walk_R2[j] += thisR2;
            walk_R2_store[i][j] = thisR2;
            // Saving coordinates
            walk_x[i][j] = indi;
            walk_y[i][j] = indj;
        }
    }
    acceptance /= (Nrun*Nrealizations);

    for(int i=0; i<Nrun; i++){
        walk_R2[i]/=Nrealizations;
        //cout << "walk_R2[" << i << "]: " << walk_R2[i] << endl;
    } // Averaging complete


    // Finding rms values
    for(int j=0; j<Nrun; j++){
        for(int i=0; i<Nrealizations; i++){
            walk_R2_rms[j] += (walk_R2_store[i][j]-walk_R2[j])*(walk_R2_store[i][j]-walk_R2[j]);
        }
        walk_R2_rms[j] = sqrt(walk_R2_rms[j]/(Nrealizations-1));
    }

    //cout << "nx: " << nx << endl << "ny: " << ny << endl;
    //cout << "Ntotblocks-1: " << Ntotblocks-1 << endl;


    //----- Partitioning -----//

    int startpoint, thispoint;
    double startx, starty, currx, curry, dx, dy, dotprod;
    for(int i=0; i<Nrealizations; i++){
        for(int j=0; j<Nsections; j++){
            // Extract starting points
            startpoint = startpoints[j];
            startx = walk_x[i][startpoint];
            starty = walk_y[i][startpoint];
            // Reset accflags:
            accflagx = 0;
            accflagy = 0;
            for(int k=1; k<len_sections; k++){ // Start at 1 since first point in array will be 0 anyways
                thispoint = startpoint+k;
                if(passing_x[i][thispoint]!=0){accflagx += passing_x[i][thispoint]*Nmat;} // Add + or - Nmat if a border is crossed
                if(passing_y[i][thispoint]!=0){accflagy += passing_y[i][thispoint]*Nmat;}
                currx   = walk_x[i][startpoint+k]+accflagx;
                curry   = walk_y[i][startpoint+k]+accflagy;
                dx      = currx-startx;
                dy      = curry-starty;
                dotprod = dx*dx+dy*dy;
                walk_R2_sections[j][k]          += dotprod;
                walk_R2_sections_store[i][j][k]  = dotprod;
            } // End loop over steps (len_sections)
        } // End loop over sections
    } // End loop over realizations

   // Averaging
    for(int i = 0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            walk_R2_sections[i][j] /= Nrealizations;
        } // End loop over steps (len_sections)
    } // End loop over sections

    // rms
    for(int i=0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            for(int k=0; k<Nrealizations; k++){
                walk_R2_sections_rms[i][j] = (walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j])*(walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j]);
            } // End loop over realizations
            walk_R2_sections_rms[i][j] = sqrt(walk_R2_sections_rms[i][j]/(Nrealizations-1));
        } // End loop over steps (len_sections)
    } // End loop over sections


    //----- Writing to file -----//
    //string outfilenamePrefix;


    ofstream outFile;
    char *filename = new char[100000]; // Possibly a long file name
    sprintf(filename, "PBC/sigma%i_d%i/Nblocks%i/pot_sigma%.3f_exp%.3f_factor%.3f_beta%.1f_R2_Nsteps%i_Nreal%i_printevery%i", intsigma, intd, Nblocks, sigma, power, soften, beta, Nrun, Nrealizations, printevery);
    //sprintf(filename, "Nsteps%i_Nreal%i_Nblocks%i_intsigma%i_intd%i_pot_exp%f_sigma%f_factor%f_PBC_R2_basic_mc", Nrun, Nrealizations, Nblocks, intsigma, intd, power, sigma, soften);
    outFile.open(filename);
    delete filename;


    for(int i=0; i<Nrun; i++){
        if(i%printevery==0){
           outFile << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl;
        }
        //outFile << std::setprecision(std::numeric_limits<double>::digits10+1) << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl; // Is this long printing really neccessary? // Possibly (probably?) gonna use this later.
    }
    outFile.close();

    ofstream poutFile;
    char *pfilename = new char[100000]; // Possibly a long file name
    sprintf(pfilename, "PBC/sigma%i_d%i/Nblocks%i/pot_sigma%.3f_exp%.3f_factor%.3f_beta%.1f_R2_Nsteps%i_Nreal%i_Npart%i_printevery%i", intsigma, intd, Nblocks, sigma, power, soften, beta, Nrun, Nrealizations, Nsections, printevery);
    //sprintf(pfilename, "Nsteps%i_Nreal%i_Nblocks%i_intsigma%i_intd%i_Npart%i_pot_exp%f_sigma%f_factor%f_PBC_R2_basic_mc", Nrun, Nrealizations, Nblocks, intsigma, intd, Nsections, power, sigma, soften);
    poutFile.open(pfilename);
    delete pfilename;

    for(int i=0; i<Nsections; i++){
        poutFile << "--- Section " << i <<" ---" << endl;
        for(int j=0; j<len_sections; j++){
            if(j%printevery==0){
                poutFile << j << " " << walk_R2_sections[i][j] << " " << walk_R2_sections_rms[i][j] << endl;
            }
        } // End loop over steps (len_sections)
    } // End loop over sections
    poutFile.close();

    if(pbc_codetesting){
        ofstream coordFile;
        char *cfilename = new char[100000]; // Possibly a long file name
        sprintf(cfilename, "PBC/sigma%i_d%i/Nblocks%i/pot_sigma%.3f_exp%.3f_factor%.3f_beta%.1f_R2_Nsteps%i_Nreal%i_coords", intsigma, intd, Nblocks, sigma, power, soften, beta, Nrun, Nrealizations);
        //sprintf(cfilename, "Nsteps%i_Nreal%i_Npart%i_pot_exp%f_sigma%f_factor%f_PBC_R2_basic_mc_coords", Nrun, Nrealizations, Nsections, power, sigma, soften);
        coordFile.open(cfilename);
        delete cfilename;

        for(int i=0; i<Nrealizations; i++){
            coordFile << "-- Realization " << i << "--" << endl;
            for(int j=0; j<Nrun; j++){
                coordFile << j << " " << walk_x[i][j] << " " << walk_y[i][j] << " ; " << passing_x[i][j] << " " << passing_y[i][j] << endl;
            }
        }
        coordFile.close();
    }


    cout << "acceptance rate: " << acceptance << endl;
    cout << "consecutive_locked: " << consecutive_locked << endl;
    cout << "e_curr: " << e_curr << endl;
}

void nearwall_rw(int maxdistance, int Nrun, int Nrealizations, int Nsections, int printevery)
{
    //---- Partitioning set up ----//
    vector<int> startpoints = give_startpoints(Nrun, Nsections);
    int len_sections = Nrun/Nsections; // I have already done this. Maybe having the function is not that handy...

    //---- Set up of run ----//
    // Random number distributions
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution_x(-50, 50); // We can be almost anywhere in the x-direction
    std::uniform_int_distribution<int> distribution_y(0, maxdistance);
    std::uniform_int_distribution<int> steps_distr(0,1); // Move: +/-1. Should I be able to keep still too?
    std::uniform_int_distribution<int> dir_distr(0,1);

    // Variables
    bool moveit;
    //double thisR2;    // Is this a conflict? Should this be an int?
    int starti, startj, indi, indj, nextindi, nextindj, dir, step, thisR2, nx, ny;
    vector<double> walk_R2(Nrun);                                           // Average R^2
    vector<double> walk_R2_rms(Nrun);                                       // RMS R^2
    vector<vector<double>> walk_R2_store(Nrealizations, vector<double>(Nrun)); // For calculation of the standard deviation
    vector<vector<double>> walk_x(Nrealizations, vector<double>(Nrun));        // Useful in case I want to look at the walks for different starting points (post-processing)
    vector<vector<double>> walk_y(Nrealizations, vector<double>(Nrun));
    vector<vector<double>> walk_R2_sections(Nsections, vector<double>(len_sections));                                                      // Sections
    vector<vector<double>> walk_R2_sections_rms(Nsections, vector<double>(len_sections));                                                  // Sections, rms
    vector<vector<vector<double>>> walk_R2_sections_store(Nrealizations, vector<vector<double>>(Nsections, vector<double>(len_sections))); // For calculation of the standard deviation


    //---- Performing walk ----//
    nx = 0; ny = 0;
    for(int i=0; i<Nrealizations; i++){
        // Start the walk by placing the bead in a random position
        // and check that the position is not occupied
        indi = distribution_x(generator); // Random numbers this way (from <random>) is a bit better than rand() from cstdlib, it seems.
        indj = distribution_y(generator);
        starti = indi;
        startj = indj;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl;
        for(int j=0; j<Nrun; j++){
            moveit = false;
            while(!moveit){    // Should I let the walker stand still if the move is not accepted?
                step = 2*steps_distr(generator)-1;
                //cout << "step: " << step << endl;
                dir = dir_distr(generator);
                //cout << "dir: " << dir << endl;
                if(dir==0){
                    //cout << "dir x chosen" << endl;
                    nextindi = indi + step;
                    nextindj = indj;
                    moveit=true;
                        nx++;
                }
                else{
                    //cout << "dir y chosen" << endl;
                    nextindi = indi;
                    nextindj = indj + step;
                    if(nextindj>-1){ // We don't hit the wall
                        moveit=true;
                        ny++;
                    }
                }
            }
            // Move the walker
            indi = nextindi;
            indj = nextindj;
            // store data somewhere
            // Average moves
            thisR2 = (indi-starti)*(indi-starti) + (indj-startj)*(indj-startj);
            walk_R2[j] += thisR2;
            walk_R2_store[i][j] = thisR2;
            // Saving coordinates
            walk_x[i][j] = indi;
            walk_y[i][j] = indj;
        }
        //cout << "indi: " << indi << ", indj: " << indj << endl;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl << "-------------------------------------" << endl;
    }

    for(int i=0; i<Nrun; i++){
        walk_R2[i]/=Nrealizations;
        //cout << "walk_R2[" << i << "]: " << walk_R2[i] << endl;
    } // Averaging complete


    // Finding rms values
    for(int j=0; j<Nrun; j++){
        for(int i=0; i<Nrealizations; i++){
            walk_R2_rms[j] += (walk_R2_store[i][j]-walk_R2[j])*(walk_R2_store[i][j]-walk_R2[j]);
        }
        walk_R2_rms[j] = sqrt(walk_R2_rms[j]/(Nrealizations-1));
    }

    //cout << "nx: " << nx << endl << "ny: " << ny << endl;
    //cout << "Nmat-1: " << Nmat-1 << endl;


    //----- Partitioning -----//

    int startpoint;
    double startx, starty, currx, curry, dx, dy, dotprod;
    for(int i=0; i<Nrealizations; i++){
        for(int j=0; j<Nsections; j++){
            // Extract starting points
            startpoint = startpoints[j];
            startx = walk_x[i][startpoint];
            starty = walk_y[i][startpoint];
            for(int k=1; k<len_sections; k++){ // Start at 1 since first point in array will be 0 anyways
                currx   = walk_x[i][startpoint+k];
                curry   = walk_y[i][startpoint+k];
                dx      = currx-startx;
                dy      = curry-starty;
                dotprod = dx*dx+dy*dy;
                walk_R2_sections[j][k]          += dotprod;
                walk_R2_sections_store[i][j][k]  = dotprod;
            } // End loop over steps (len_sections)
        } // End loop over sections
    } // End loop over realizations

   // Averaging
    for(int i = 0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            walk_R2_sections[i][j] /= Nrealizations;
        } // End loop over steps (len_sections)
    } // End loop over sections

    // rms
    for(int i=0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            for(int k=0; k<Nrealizations; k++){
                walk_R2_sections_rms[i][j] = (walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j])*(walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j]);
            } // End loop over realizations
            walk_R2_sections_rms[i][j] = sqrt(walk_R2_sections_rms[i][j]/(Nrealizations-1));
        } // End loop over steps (len_sections)
    } // End loop over sections


    //----- Writing to file -----//
    //string outfilenamePrefix;

    ofstream outFile;
    char *filename = new char[100000]; // Possibly a long file name
    sprintf(filename, "nearwallrw_R2_Nsteps%i_Nreal%i_maxstartdist%i_printevery%i", Nrun, Nrealizations, maxdistance, printevery);
    outFile.open(filename);
    delete filename;


    for(int i=0; i<Nrun; i++){
        if(i%printevery==0){
            outFile << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl;
        }
    }
    outFile.close();

    ofstream poutFile;
    char *pfilename = new char[100000]; // Possibly a long file name
    sprintf(pfilename, "nearwallrw_R2_Nsteps%i_Nreal%i_maxstartdist%i_Npart%i_printevery%i", Nrun, Nrealizations, maxdistance, Nsections, printevery);
    poutFile.open(pfilename);
    delete pfilename;

    for(int i=0; i<Nsections; i++){
        poutFile << "--- Section " << i <<" ---" << endl;
        for(int j=0; j<len_sections; j++){
            if(j%printevery==0){
                poutFile << j << " " << walk_R2_sections[i][j] << " " << walk_R2_sections_rms[i][j] << endl;
            }
        } // End loop over steps (len_sections)
    } // End loop over sections
    poutFile.close();
}

void nearwall_mc(int maxdistance, int Nrun, int Nrealizations, int Nsections, int printevery)
{
    //---- Partitioning set up ----//
    vector<int> startpoints = give_startpoints(Nrun, Nsections);
    int len_sections = Nrun/Nsections; // I have already done this. Maybe having the function is not that handy...

    //---- Set up of run ----//
    // Random number distributions
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution_x(-50, 50); // We can be almost anywhere in the x-direction
    std::uniform_int_distribution<int> distribution_y(0, maxdistance);
    std::uniform_int_distribution<int> steps_distr(0,1); // Move: +/-1. Should I be able to keep still too?
    std::uniform_int_distribution<int> dir_distr(0,1);

    // Variables
    bool moveit;
    //double thisR2;    // Is this a conflict? Should this be an int?
    int starti, startj, indi, indj, nextindi, nextindj, dir, step, thisR2, nx, ny;
    vector<double> walk_R2(Nrun);                                           // Average R^2
    vector<double> walk_R2_rms(Nrun);                                       // RMS R^2
    vector<vector<double>> walk_R2_store(Nrealizations, vector<double>(Nrun)); // For calculation of the standard deviation
    vector<vector<double>> walk_x(Nrealizations, vector<double>(Nrun));        // Useful in case I want to look at the walks for different starting points (post-processing)
    vector<vector<double>> walk_y(Nrealizations, vector<double>(Nrun));
    vector<vector<double>> walk_R2_sections(Nsections, vector<double>(len_sections));                                                      // Sections
    vector<vector<double>> walk_R2_sections_rms(Nsections, vector<double>(len_sections));                                                  // Sections, rms
    vector<vector<vector<double>>> walk_R2_sections_store(Nrealizations, vector<vector<double>>(Nsections, vector<double>(len_sections))); // For calculation of the standard deviation


    //---- Performing walk ----//
    nx = 0; ny = 0;
    for(int i=0; i<Nrealizations; i++){
        // Start the walk by placing the bead in a random position
        // and check that the position is not occupied
        indi = distribution_x(generator); // Random numbers this way (from <random>) is a bit better than rand() from cstdlib, it seems.
        indj = distribution_y(generator);
        starti = indi;
        startj = indj;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl;
        for(int j=0; j<Nrun; j++){
            step = 2*steps_distr(generator)-1;
            //cout << "step: " << step << endl;
            dir = dir_distr(generator);
            //cout << "dir: " << dir << endl;
            if(dir==0){
                //cout << "dir x chosen" << endl;
                nextindi = indi + step;
                nextindj = indj;
                moveit=true;
                nx++;
            }
            else{
                //cout << "dir y chosen" << endl;
                nextindi = indi;
                nextindj = indj + step;
                if(nextindj>-1){ // We don't hit the wall
                    ny++;
                }
                else{
                    nextindj = indj;
                }
            }
            // Move the walker
            indi = nextindi;
            indj = nextindj;
            // store data somewhere
            // Average moves
            thisR2 = (indi-starti)*(indi-starti) + (indj-startj)*(indj-startj);
            walk_R2[j] += thisR2;
            walk_R2_store[i][j] = thisR2;
            // Saving coordinates
            walk_x[i][j] = indi;
            walk_y[i][j] = indj;
        }
        //cout << "indi: " << indi << ", indj: " << indj << endl;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl << "-------------------------------------" << endl;
    }

    for(int i=0; i<Nrun; i++){
        walk_R2[i]/=Nrealizations;
        //cout << "walk_R2[" << i << "]: " << walk_R2[i] << endl;
    } // Averaging complete


    // Finding rms values
    for(int j=0; j<Nrun; j++){
        for(int i=0; i<Nrealizations; i++){
            walk_R2_rms[j] += (walk_R2_store[i][j]-walk_R2[j])*(walk_R2_store[i][j]-walk_R2[j]);
        }
        walk_R2_rms[j] = sqrt(walk_R2_rms[j]/(Nrealizations-1));
    }

    //cout << "nx: " << nx << endl << "ny: " << ny << endl;
    //cout << "Nmat-1: " << Nmat-1 << endl;


    //----- Partitioning -----//

    int startpoint;
    double startx, starty, currx, curry, dx, dy, dotprod;
    for(int i=0; i<Nrealizations; i++){
        for(int j=0; j<Nsections; j++){
            // Extract starting points
            startpoint = startpoints[j];
            startx = walk_x[i][startpoint];
            starty = walk_y[i][startpoint];
            for(int k=1; k<len_sections; k++){ // Start at 1 since first point in array will be 0 anyways
                currx   = walk_x[i][startpoint+k];
                curry   = walk_y[i][startpoint+k];
                dx      = currx-startx;
                dy      = curry-starty;
                dotprod = dx*dx+dy*dy;
                walk_R2_sections[j][k]          += dotprod;
                walk_R2_sections_store[i][j][k]  = dotprod;
            } // End loop over steps (len_sections)
        } // End loop over sections
    } // End loop over realizations

   // Averaging
    for(int i = 0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            walk_R2_sections[i][j] /= Nrealizations;
        } // End loop over steps (len_sections)
    } // End loop over sections

    // rms
    for(int i=0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            for(int k=0; k<Nrealizations; k++){
                walk_R2_sections_rms[i][j] = (walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j])*(walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j]);
            } // End loop over realizations
            walk_R2_sections_rms[i][j] = sqrt(walk_R2_sections_rms[i][j]/(Nrealizations-1));
        } // End loop over steps (len_sections)
    } // End loop over sections


    //----- Writing to file -----//
    //string outfilenamePrefix;

    ofstream outFile;
    char *filename = new char[100000]; // Possibly a long file name
    sprintf(filename, "nearwallmc_R2_Nsteps%i_Nreal%i_maxstartdist%i_printevery%i", Nrun, Nrealizations, maxdistance, printevery);
    outFile.open(filename);
    delete filename;


    for(int i=0; i<Nrun; i++){
        if(i%printevery==0){
            outFile << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl;
        }
    }
    outFile.close();

    ofstream poutFile;
    char *pfilename = new char[100000]; // Possibly a long file name
    sprintf(pfilename, "nearwallmc_R2_Nsteps%i_Nreal%i_maxstartdist%i_Npart%i_printevery%i", Nrun, Nrealizations, maxdistance, Nsections, printevery);
    poutFile.open(pfilename);
    delete pfilename;

    for(int i=0; i<Nsections; i++){
        poutFile << "--- Section " << i <<" ---" << endl;
        for(int j=0; j<len_sections; j++){
            if(j%printevery==0){
                poutFile << j << " " << walk_R2_sections[i][j] << " " << walk_R2_sections_rms[i][j] << endl;
            }
        } // End loop over steps (len_sections)
    } // End loop over sections
    poutFile.close();
}


void randomwalk_1D(int Nrun, int Nrealizations, int Nsections, int printevery)
{
    //---- Partitioning set up ----//
    vector<int> startpoints = give_startpoints(Nrun, Nsections);
    int len_sections = Nrun/Nsections; // I have already done this. Maybe having the function is not that handy...

    //---- Set up of run ----//
    // Random number distributions
    std::default_random_engine generator;
    std::uniform_int_distribution<int> steps_distr(0,1); // Move: +/-1.

    // Variables
    bool moveit;
    //double thisR2;    // Is this a conflict? Should this be an int?
    int x, step, thisR2;
    vector<double> walk_R2(Nrun);                                           // Average R^2
    vector<double> walk_R2_rms(Nrun);                                       // RMS R^2
    vector<vector<double>> walk_R2_store(Nrealizations, vector<double>(Nrun)); // For calculation of the standard deviation
    vector<vector<double>> walk_x(Nrealizations, vector<double>(Nrun));        // Useful in case I want to look at the walks for different starting points (post-processing)
    vector<vector<double>> walk_R2_sections(Nsections, vector<double>(len_sections));                                                      // Sections
    vector<vector<double>> walk_R2_sections_rms(Nsections, vector<double>(len_sections));                                                  // Sections, rms
    vector<vector<vector<double>>> walk_R2_sections_store(Nrealizations, vector<vector<double>>(Nsections, vector<double>(len_sections))); // For calculation of the standard deviation


    //---- Performing walk ----//
    for(int i=0; i<Nrealizations; i++){
        // Start the walk by placing the bead in a random position
        // and check that the position is not occupied
        x = 0;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl;
        for(int j=0; j<Nrun; j++){
            // Draw: Move forwards or backwards
            step = 2*steps_distr(generator)-1;
            // Move the walker
            x += step;
            // store data somewhere
            // Average moves
            thisR2 = x*x;
            walk_R2[j] += thisR2;
            walk_R2_store[i][j] = thisR2;
            // Saving coordinates
            walk_x[i][j] = x;
        }
        //cout << "indi: " << indi << ", indj: " << indj << endl;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl << "-------------------------------------" << endl;
    }

    for(int i=0; i<Nrun; i++){
        walk_R2[i]/=Nrealizations;
        //cout << "walk_R2[" << i << "]: " << walk_R2[i] << endl;
    } // Averaging complete


    // Finding rms values
    for(int j=0; j<Nrun; j++){
        for(int i=0; i<Nrealizations; i++){
            walk_R2_rms[j] += (walk_R2_store[i][j]-walk_R2[j])*(walk_R2_store[i][j]-walk_R2[j]);
        }
        walk_R2_rms[j] = sqrt(walk_R2_rms[j]/(Nrealizations-1));
    }

    //cout << "nx: " << nx << endl << "ny: " << ny << endl;
    //cout << "Nmat-1: " << Nmat-1 << endl;


    //----- Partitioning -----//

    int startpoint;
    double startx, currx, dx, dotprod;
    for(int i=0; i<Nrealizations; i++){
        for(int j=0; j<Nsections; j++){
            // Extract starting points
            startpoint = startpoints[j];
            startx = walk_x[i][startpoint];
            for(int k=1; k<len_sections; k++){ // Start at 1 since first point in array will be 0 anyways
                currx   = walk_x[i][startpoint+k];
                dx      = currx-startx;
                dotprod = dx*dx;
                walk_R2_sections[j][k]          += dotprod;
                walk_R2_sections_store[i][j][k]  = dotprod;
            } // End loop over steps (len_sections)
        } // End loop over sections
    } // End loop over realizations

   // Averaging
    for(int i = 0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            walk_R2_sections[i][j] /= Nrealizations;
        } // End loop over steps (len_sections)
    } // End loop over sections

    // rms
    for(int i=0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            for(int k=0; k<Nrealizations; k++){
                walk_R2_sections_rms[i][j] = (walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j])*(walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j]);
            } // End loop over realizations
            walk_R2_sections_rms[i][j] = sqrt(walk_R2_sections_rms[i][j]/(Nrealizations-1));
        } // End loop over steps (len_sections)
    } // End loop over sections


    //----- Writing to file -----//
    //string outfilenamePrefix;

    ofstream outFile;
    char *filename = new char[100000]; // Possibly a long file name
    sprintf(filename, "rw1D_R2_Nsteps%i_Nreal%i_printevery%i", Nrun, Nrealizations, printevery);
    outFile.open(filename);
    delete filename;


    for(int i=0; i<Nrun; i++){
        if(i%printevery==0){
            outFile << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl;
        }
    }
    outFile.close();

    ofstream poutFile;
    char *pfilename = new char[100000]; // Possibly a long file name
    sprintf(pfilename, "rw1D_R2_Nsteps%i_Nreal%i_Npart%i_printevery%i", Nrun, Nrealizations, Nsections, printevery);
    poutFile.open(pfilename);
    delete pfilename;

    for(int i=0; i<Nsections; i++){
        poutFile << "--- Section " << i <<" ---" << endl;
        for(int j=0; j<len_sections; j++){
            if(j%printevery==0){
                poutFile << j << " " << walk_R2_sections[i][j] << " " << walk_R2_sections_rms[i][j] << endl;
            }
        } // End loop over steps (len_sections)
    } // End loop over sections
    poutFile.close();
}


//---------------- Not so flexible functions ----------------//
void matrixwalk_hard(int Nblocks, int blocksize, int Nrun, int Nrealizations, int Nsections)
{
    //---- Matrix set up ----//
    int Ntotblocks, Nmat;
    Nmat = Nblocks*blocksize; // Empty and filled space is of the same size. Maybe a temporary solution?
    Ntotblocks = (2*Nblocks+1)*blocksize;
    vector<vector<int>> walkmat(Ntotblocks, vector<int>(Ntotblocks));

    // Setting matrix blocks
    // Must be a handier way to do this?:
    for(int i=0; i<Nblocks; i++){
        for(int j=0; j<Nblocks; j++){
            for(int k=0; k<blocksize; k++){
                for(int l=0; l<blocksize; l++){
                    walkmat[(2*i+1)*blocksize+k][(2*j+1)*blocksize+l] = 1;
                }
            }
        }
    }

    //---- Partitioning set up ----//
    vector<int> startpoints = give_startpoints(Nrun, Nsections);
    int len_sections = Nrun/Nsections; // I have already done this. Maybe having the function is not that handy...

    //---- Set up of run ----//
    // Random number distributions
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,Ntotblocks-1);
    std::uniform_int_distribution<int> steps_distr(0,1); // Move: +/-1. Should I be able to keep still too?
    std::uniform_int_distribution<int> dir_distr(0,1);

    // Variables
    bool free, moveit;
    //double thisR2;    // Is this a conflict? Should this be an int?
    int starti, startj, indi, indj, nextindi, nextindj, dir, step, thisR2, nx, ny;
    vector<double> walk_R2(Nrun);                                           // Average R^2
    vector<double> walk_R2_rms(Nrun);                                       // RMS R^2
    vector<vector<double>> walk_R2_store(Nrealizations, vector<double>(Nrun)); // For calculation of the standard deviation
    vector<vector<double>> walk_x(Nrealizations, vector<double>(Nrun));        // Useful in case I want to look at the walks for different starting points (post-processing)
    vector<vector<double>> walk_y(Nrealizations, vector<double>(Nrun));
    vector<vector<double>> walk_R2_sections(Nsections, vector<double>(len_sections));                                                      // Sections
    vector<vector<double>> walk_R2_sections_rms(Nsections, vector<double>(len_sections));                                                  // Sections, rms
    vector<vector<vector<double>>> walk_R2_sections_store(Nrealizations, vector<vector<double>>(Nsections, vector<double>(len_sections))); // For calculation of the standard deviation


    //---- Performing walk ----//
    nx = 0; ny = 0;
    for(int i=0; i<Nrealizations; i++){
        // Start the walk by placing the bead in a random position
        // and check that the position is not occupied
        free = false;
        while(!free){
            indi = distribution(generator); //location();//indi = rand() % Nmat; // Random numbers this way (from <random>) is a bit better than rand() from cstdlib, it seems.
            indj = distribution(generator); //location();//indj = rand() % Nmat;
            //cout << "walkmat["<<indi<<"]["<<indj<< "]: " << walkmat[indi][indj] << endl;
            if(walkmat[indi][indj]==0){free=true;} // Only place the bead on a free site
        }
        starti = indi;
        startj = indj;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl;
        for(int j=0; j<Nrun; j++){
            moveit = false;
            while(!moveit){    // Should I let the walker stand still if the move is not accepted?
                step = 2*steps_distr(generator)-1;
                //cout << "step: " << step << endl;
                dir = dir_distr(generator);
                //cout << "dir: " << dir << endl;
                if(dir==0){
                    //cout << "dir x chosen" << endl;
                    nextindi = indi + step;
                    nextindj = indj;
                    if(nextindi<Ntotblocks && nextindi>-1){
                        if(walkmat[nextindi][nextindj]==0){
                            moveit=true;
                            nx++;
                        }
                    }
                }
                else{
                    //cout << "dir y chosen" << endl;
                    nextindi = indi;
                    nextindj = indj + step;
                    if(nextindj<Ntotblocks && nextindj>-1){ // If the site exists
                        if(walkmat[nextindi][nextindj]==0){
                            moveit=true;
                            ny++;
                        }
                    }
                }
            }
            // Move the walker
            indi = nextindi;
            indj = nextindj;
            // store data somewhere
            // Average moves
            thisR2 = (indi-starti)*(indi-starti) + (indj-startj)*(indj-startj);
            walk_R2[j] += thisR2;
            walk_R2_store[i][j] = thisR2;
            // Saving coordinates
            walk_x[i][j] = indi;
            walk_y[i][j] = indj;
        }
        //cout << "indi: " << indi << ", indj: " << indj << endl;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl << "-------------------------------------" << endl;
    }

    for(int i=0; i<Nrun; i++){
        walk_R2[i]/=Nrealizations;
        //cout << "walk_R2[" << i << "]: " << walk_R2[i] << endl;
    } // Averaging complete


    // Finding rms values
    for(int j=0; j<Nrun; j++){
        for(int i=0; i<Nrealizations; i++){
            walk_R2_rms[j] += (walk_R2_store[i][j]-walk_R2[j])*(walk_R2_store[i][j]-walk_R2[j]);
        }
        walk_R2_rms[j] = sqrt(walk_R2_rms[j]/(Nrealizations-1));
    }

    //cout << "nx: " << nx << endl << "ny: " << ny << endl;
    //cout << "Ntotblocks-1: " << Ntotblocks-1 << endl;


    //----- Partitioning -----//

    int startpoint;
    double startx, starty, currx, curry, dx, dy, dotprod;
    for(int i=0; i<Nrealizations; i++){
        for(int j=0; j<Nsections; j++){
            // Extract starting points
            startpoint = startpoints[j];
            startx = walk_x[i][startpoint];
            starty = walk_y[i][startpoint];
            for(int k=1; k<len_sections; k++){ // Start at 1 since first point in array will be 0 anyways
                currx   = walk_x[i][startpoint+k];
                curry   = walk_y[i][startpoint+k];
                dx      = currx-startx;
                dy      = curry-starty;
                dotprod = dx*dx+dy*dy;
                walk_R2_sections[j][k]          += dotprod;
                walk_R2_sections_store[i][j][k]  = dotprod;
            } // End loop over steps (len_sections)
        } // End loop over sections
    } // End loop over realizations

   // Averaging
    for(int i = 0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            walk_R2_sections[i][j] /= Nrealizations;
        } // End loop over steps (len_sections)
    } // End loop over sections

    // rms
    for(int i=0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            for(int k=0; k<Nrealizations; k++){
                walk_R2_sections_rms[i][j] = (walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j])*(walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j]);
            } // End loop over realizations
            walk_R2_sections_rms[i][j] = sqrt(walk_R2_sections_rms[i][j]/(Nrealizations-1));
        } // End loop over steps (len_sections)
    } // End loop over sections


    //----- Writing to file -----//
    //string outfilenamePrefix;


    ofstream outFile;
    char *filename = new char[100000]; // Possibly a long file name
    sprintf(filename, "Nsteps%i_Nreal%i_hardpot_R2_basic", Nrun, Nrealizations);
    outFile.open(filename);
    delete filename;


    for(int i=0; i<Nrun; i++){
        outFile << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl;
        //outFile << std::setprecision(std::numeric_limits<double>::digits10+1) << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl; // Is this long printing really neccessary? // Possibly (probably?) gonna use this later.
    }
    outFile.close();

    ofstream poutFile;
    char *pfilename = new char[100000]; // Possibly a long file name
    sprintf(pfilename, "Nsteps%i_Nreal%i_Npart%i_hardpot_R2_basic", Nrun, Nrealizations, Nsections);
    poutFile.open(pfilename);
    delete pfilename;

    for(int i=0; i<Nsections; i++){
        poutFile << "--- Section " << i <<" ---" << endl;
        for(int j=0; j<len_sections; j++){
            poutFile << j << " " << walk_R2_sections[i][j] << " " << walk_R2_sections_rms[i][j] << endl;
        } // End loop over steps (len_sections)
    } // End loop over sections
    poutFile.close();
}

void matrixmc_hard(int Nblocks, int blocksize, int Nrun, int Nrealizations, int Nsections)
{
    //---- Matrix set up ----//
    int Ntotblocks, Nmat;
    Nmat = Nblocks*blocksize; // Empty and filled space is of the same size. Maybe a temporary solution?
    Ntotblocks = (2*Nblocks+1)*blocksize;
    vector<vector<int>> walkmat(Ntotblocks, vector<int>(Ntotblocks));

    // Setting matrix blocks
    // Must be a handier way to do this?:
    for(int i=0; i<Nblocks; i++){
        for(int j=0; j<Nblocks; j++){
            for(int k=0; k<blocksize; k++){
                for(int l=0; l<blocksize; l++){
                    walkmat[(2*i+1)*blocksize+k][(2*j+1)*blocksize+l] = 1;
                }
            }
        }
    }

    //---- Partitioning set up ----//
    vector<int> startpoints = give_startpoints(Nrun, Nsections);
    int len_sections = Nrun/Nsections; // I have already done this. Maybe having the function is not that handy...

    //---- Set up of run ----//
    // Random number distributions
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,Ntotblocks-1);
    std::uniform_int_distribution<int> steps_distr(0,1); // Move: +/-1. Should I be able to keep still too?
    std::uniform_int_distribution<int> dir_distr(0,1);

    // Variables
    bool free;
    //double thisR2;    // Is this a conflict? Should this be an int?
    int starti, startj, indi, indj, nextindi, nextindj, dir, step, thisR2, nx, ny;
    vector<double> walk_R2(Nrun);                                           // Average R^2
    vector<double> walk_R2_rms(Nrun);                                       // RMS R^2
    vector<vector<double>> walk_R2_store(Nrealizations, vector<double>(Nrun)); // For calculation of the standard deviation
    vector<vector<double>> walk_x(Nrealizations, vector<double>(Nrun));        // Useful in case I want to look at the walks for different starting points (post-processing)
    vector<vector<double>> walk_y(Nrealizations, vector<double>(Nrun));
    vector<vector<double>> walk_R2_sections(Nsections, vector<double>(len_sections));                                                      // Sections
    vector<vector<double>> walk_R2_sections_rms(Nsections, vector<double>(len_sections));                                                  // Sections, rms
    vector<vector<vector<double>>> walk_R2_sections_store(Nrealizations, vector<vector<double>>(Nsections, vector<double>(len_sections))); // For calculation of the standard deviation


    //---- Performing walk ----//
    nx = 0; ny = 0;
    for(int i=0; i<Nrealizations; i++){
        // Start the walk by placing the bead in a random position
        // and check that the position is not occupied
        free = false;
        while(!free){
            indi = distribution(generator); //location();//indi = rand() % Nmat; // Random numbers this way (from <random>) is a bit better than rand() from cstdlib, it seems.
            indj = distribution(generator); //location();//indj = rand() % Nmat;
            //cout << "walkmat["<<indi<<"]["<<indj<< "]: " << walkmat[indi][indj] << endl;
            if(walkmat[indi][indj]==0){free=true;} // Only place the bead on a free site
        }
        starti = indi;
        startj = indj;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl;
        for(int j=0; j<Nrun; j++){
            step = 2*steps_distr(generator)-1;
            //cout << "step: " << step << endl;
            dir = dir_distr(generator);
            //cout << "dir: " << dir << endl;
            if(dir==0){
                //cout << "dir x chosen" << endl;
                nextindi = indi + step;
                nextindj = indj;
                if(nextindi<Ntotblocks && nextindi>-1){
                    if(walkmat[nextindi][nextindj]==0){
                        nx++;
                    }
                    else{
                        nextindi = indi;
                        nextindj = indj;
                    }
                }
                else{
                    nextindi = indi;
                    nextindj = indj;
                }
            }
            else{
                //cout << "dir y chosen" << endl;
                nextindi = indi;
                nextindj = indj + step;
                if(nextindj<Ntotblocks && nextindj>-1){ // If the site exists
                    if(walkmat[nextindi][nextindj]==0){
                        ny++;
                    }
                    else{
                        nextindi = indi;
                        nextindj = indj;
                    }
                }
                else{
                    nextindi = indi;
                    nextindj = indj;
                }
            }
            // Move the walker
            indi = nextindi;
            indj = nextindj;
            // store data somewhere
            // Average moves
            thisR2 = (indi-starti)*(indi-starti) + (indj-startj)*(indj-startj);
            walk_R2[j] += thisR2;
            walk_R2_store[i][j] = thisR2;
            // Saving coordinates
            walk_x[i][j] = indi;
            walk_y[i][j] = indj;
        }
    }

    for(int i=0; i<Nrun; i++){
        walk_R2[i]/=Nrealizations;
        //cout << "walk_R2[" << i << "]: " << walk_R2[i] << endl;
    } // Averaging complete


    // Finding rms values
    for(int j=0; j<Nrun; j++){
        for(int i=0; i<Nrealizations; i++){
            walk_R2_rms[j] += (walk_R2_store[i][j]-walk_R2[j])*(walk_R2_store[i][j]-walk_R2[j]);
        }
        walk_R2_rms[j] = sqrt(walk_R2_rms[j]/(Nrealizations-1));
    }

    //cout << "nx: " << nx << endl << "ny: " << ny << endl;
    //cout << "Ntotblocks-1: " << Ntotblocks-1 << endl;


    //----- Partitioning -----//

    int startpoint;
    double startx, starty, currx, curry, dx, dy, dotprod;
    for(int i=0; i<Nrealizations; i++){
        for(int j=0; j<Nsections; j++){
            // Extract starting points
            startpoint = startpoints[j];
            startx = walk_x[i][startpoint];
            starty = walk_y[i][startpoint];
            for(int k=1; k<len_sections; k++){ // Start at 1 since first point in array will be 0 anyways
                currx   = walk_x[i][startpoint+k];
                curry   = walk_y[i][startpoint+k];
                dx      = currx-startx;
                dy      = curry-starty;
                dotprod = dx*dx+dy*dy;
                walk_R2_sections[j][k]          += dotprod;
                walk_R2_sections_store[i][j][k]  = dotprod;
            } // End loop over steps (len_sections)
        } // End loop over sections
    } // End loop over realizations

   // Averaging
    for(int i = 0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            walk_R2_sections[i][j] /= Nrealizations;
        } // End loop over steps (len_sections)
    } // End loop over sections

    // rms
    for(int i=0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            for(int k=0; k<Nrealizations; k++){
                walk_R2_sections_rms[i][j] = (walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j])*(walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j]);
            } // End loop over realizations
            walk_R2_sections_rms[i][j] = sqrt(walk_R2_sections_rms[i][j]/(Nrealizations-1));
        } // End loop over steps (len_sections)
    } // End loop over sections


    //----- Writing to file -----//
    //string outfilenamePrefix;


    ofstream outFile;
    char *filename = new char[100000]; // Possibly a long file name
    sprintf(filename, "Nsteps%i_Nreal%i_hardpot_R2_basic_mc", Nrun, Nrealizations);
    outFile.open(filename);
    delete filename;


    for(int i=0; i<Nrun; i++){
        outFile << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl;
        //outFile << std::setprecision(std::numeric_limits<double>::digits10+1) << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl; // Is this long printing really neccessary? // Possibly (probably?) gonna use this later.
    }
    outFile.close();

    ofstream poutFile;
    char *pfilename = new char[100000]; // Possibly a long file name
    sprintf(pfilename, "Nsteps%i_Nreal%i_Npart%i_hardpot_R2_basic_mc", Nrun, Nrealizations, Nsections);
    poutFile.open(pfilename);
    delete pfilename;

    for(int i=0; i<Nsections; i++){
        poutFile << "--- Section " << i <<" ---" << endl;
        for(int j=0; j<len_sections; j++){
            poutFile << j << " " << walk_R2_sections[i][j] << " " << walk_R2_sections_rms[i][j] << endl;
        } // End loop over steps (len_sections)
    } // End loop over sections
    poutFile.close();
}


void matrixmc_potential(int Nblocks, int blocksize, int Nrun, int Nrealizations, int Nsections, double beta, double sigma, double power, double soften){
    //---- Matrix set up ----//
    // (to make it the same size as the hard potential one)
    int Ntotblocks, Nmat;
    Nmat = Nblocks*blocksize; // Empty and filled space is of the same size. Maybe a temporary solution?
    Ntotblocks = (2*Nblocks+1)*blocksize;

    //---- Potential set up ----//
    // Number of 'obstacles' given by Nblocks x Nblocks
    // Centre them evenly in a grid-like manner
    int Npots = Nblocks*Nblocks;
    double obspacing = Ntotblocks/(Nblocks+1); // Can have off-grid positions for the centres. // +1: Don't want to position the obstacles at the ends of the matrix
    vector<double> centres_x = vector<double>(Npots);
    vector<double> centres_y = vector<double>(Npots);

    cout << "obspacing: " << obspacing << endl;

    ofstream brushFile;
    char *bfilename = new char[100000]; // Possibly a long file name
    sprintf(bfilename, "Nsteps%i_Nreal%i_pot_exp%f_sigma%f_factor%f_R2_basic_mc_brush", Nrun, Nrealizations, power, sigma, soften);
    brushFile.open(bfilename);
    delete bfilename;

    int counter = 0;

    for(int i=0; i<Nblocks; i++){
        for(int j=0; j<Nblocks; j++){
            centres_x[counter] = (i+1)*obspacing;
            centres_y[counter] = (j+1)*obspacing;
            brushFile << centres_x[counter] << " " << centres_y[counter] << endl;
            counter++;
        }
    }
    brushFile.close();

    //---- Partitioning set up ----//
    vector<int> startpoints = give_startpoints(Nrun, Nsections);
    int len_sections = Nrun/Nsections; // I have already done this. Maybe having the function is not that handy...

    //---- Set up of run ----//
    // Random number distributions
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,Ntotblocks-1);
    std::uniform_int_distribution<int> steps_distr(0,1); // Move: +/-1. Should I be able to keep still too?
    std::uniform_int_distribution<int> dir_distr(0,1);
    std::uniform_real_distribution<double> accept_prob(0,1); // Probability of accepting step with higher energy

    // Variables
    bool free;
    //double thisR2;    // Is this a conflict? Should this be an int?
    int starti, startj, indi, indj, nextindi, nextindj, dir, step, thisR2, nx, ny;
    vector<double> walk_R2(Nrun);                                           // Average R^2
    vector<double> walk_R2_rms(Nrun);                                       // RMS R^2
    vector<vector<double>> walk_R2_store(Nrealizations, vector<double>(Nrun)); // For calculation of the standard deviation
    vector<vector<double>> walk_x(Nrealizations, vector<double>(Nrun));        // Useful in case I want to look at the walks for different starting points (post-processing)
    vector<vector<double>> walk_y(Nrealizations, vector<double>(Nrun));
    vector<vector<double>> walk_R2_sections(Nsections, vector<double>(len_sections));                                                      // Sections
    vector<vector<double>> walk_R2_sections_rms(Nsections, vector<double>(len_sections));                                                  // Sections, rms
    vector<vector<vector<double>>> walk_R2_sections_store(Nrealizations, vector<vector<double>>(Nsections, vector<double>(len_sections))); // For calculation of the standard deviation


    //---- Performing walk ----//
    bool last_locked;
    int consecutive_locked;
    double e_curr, e_next, de, boltzmann, acceptance, illegalmovesuggested;
    nx = 0; ny = 0; acceptance = 0;
    for(int i=0; i<Nrealizations; i++){
        // Resetting some diagnostics quantities:
        last_locked = false;
        consecutive_locked = 0;
        // Start the walk by placing the bead in a random position
        // and check that the position is not occupied
        free = false;
        while(!free){
            indi = distribution(generator); //location();//indi = rand() % Nmat; // Random numbers this way (from <random>) is a bit better than rand() from cstdlib, it seems.
            indj = distribution(generator); //location();//indj = rand() % Nmat;
            //cout << "Initiation. indi: " << indi << "; indj: " << indj << endl;
            e_curr = give_energy(indi, indj, Npots, sigma, power, soften, centres_x, centres_y);
            boltzmann = exp(-beta*e_curr);
            //cout << "walkmat["<<indi<<"]["<<indj<< "]: " << walkmat[indi][indj] << endl;
            //if(walkmat[indi][indj]==0){free=true;} // Only place the bead on a free site
            //cout << "i = " << i << "; Energy: " << e_curr << "; Boltzmann factor:" << boltzmann << endl;
            if(boltzmann>0.5){free=true;} // I've just set some arbitrary threshold so that it won't get too close to a bead.
        }
        starti = indi;
        startj = indj;
        //cout << "i = " << i << ", starti = " << starti << ", startj = " << startj << endl;
        for(int j=0; j<Nrun; j++){
            step = 2*steps_distr(generator)-1;
            //cout << "step: " << step << endl;
            dir = dir_distr(generator);
            //cout << "dir: " << dir << endl;
            if(dir==0){
                //cout << "dir x chosen" << endl;
                nextindi = indi + step;
                nextindj = indj;
                if(nextindi>(Ntotblocks-1) || nextindi<0){
                    nextindi = indi;
                    nextindj = indj;
                    illegalmovesuggested += 1;
                }
            }
            else{
                //cout << "dir y chosen" << endl;
                nextindi = indi;
                nextindj = indj + step;
                if(nextindj>(Ntotblocks-1) || nextindj<0){
                    nextindi = indi;
                    nextindj = indj;
                    illegalmovesuggested += 1;
                }
            }
            e_next    = give_energy(nextindi, nextindj, Npots, sigma, power, soften, centres_x, centres_y);
            if(e_next<e_curr){
                // Accept move
                indi = nextindi;
                indj = nextindj;
                e_curr = e_next;
                acceptance += 1;
                last_locked = false;
                consecutive_locked = 0;
                //cout << "Energy changed, e_curr: " << e_curr << endl;
            }
            else{
                de        = e_next-e_curr;
                boltzmann = exp(-beta*de);
                if(accept_prob(generator)<boltzmann){
                    //cout << "e_next: " << e_next << "; e_curr: " << e_curr << endl;
                    indi = nextindi;
                    indj = nextindj;
                    e_curr = e_next;
                    acceptance += 1;
                    last_locked = false;
                    consecutive_locked = 0;
                    //cout << "Energy changed, e_curr: " << e_curr << endl;
                }
                else{
                    if(last_locked){
                        consecutive_locked++;
                        //cout << "cout: consecutive_locked:" << consecutive_locked << "Energy not changed, e_curr: " << e_curr << "; e_next: " << e_next << "; Current position: [" << indi << "," << indj << "]" << endl;
                        /*
                        if(consecutive_locked>1000){
                            cout << "Energy not changed, e_curr: " << e_curr << "; e_next: " << e_next << "; Current position: [" << indi << "," << indj << "]" << endl;
                        }
                        */
                    }
                    last_locked = true;
                }
            }
            // Average moves
            thisR2 = (indi-starti)*(indi-starti) + (indj-startj)*(indj-startj);
            walk_R2[j] += thisR2;
            walk_R2_store[i][j] = thisR2;
            // Saving coordinates
            walk_x[i][j] = indi;
            walk_y[i][j] = indj;
        }
    }
    acceptance /= (Nrun*Nrealizations);

    for(int i=0; i<Nrun; i++){
        walk_R2[i]/=Nrealizations;
        //cout << "walk_R2[" << i << "]: " << walk_R2[i] << endl;
    } // Averaging complete


    // Finding rms values
    for(int j=0; j<Nrun; j++){
        for(int i=0; i<Nrealizations; i++){
            walk_R2_rms[j] += (walk_R2_store[i][j]-walk_R2[j])*(walk_R2_store[i][j]-walk_R2[j]);
        }
        walk_R2_rms[j] = sqrt(walk_R2_rms[j]/(Nrealizations-1));
    }

    //cout << "nx: " << nx << endl << "ny: " << ny << endl;
    //cout << "Ntotblocks-1: " << Ntotblocks-1 << endl;


    //----- Partitioning -----//

    int startpoint;
    double startx, starty, currx, curry, dx, dy, dotprod;
    for(int i=0; i<Nrealizations; i++){
        for(int j=0; j<Nsections; j++){
            // Extract starting pointsillegalmovesuggestedillegalmovesuggestedillegalmovesuggestedillegalmovesuggestedillegalmovesuggested
            startpoint = startpoints[j];
            startx = walk_x[i][startpoint];
            starty = walk_y[i][startpoint];
            for(int k=1; k<len_sections; k++){ // Start at 1 since first point in array will be 0 anyways
                currx   = walk_x[i][startpoint+k];
                curry   = walk_y[i][startpoint+k];
                dx      = currx-startx;
                dy      = curry-starty;
                dotprod = dx*dx+dy*dy;
                walk_R2_sections[j][k]          += dotprod;
                walk_R2_sections_store[i][j][k]  = dotprod;
            } // End loop over steps (len_sections)
        } // End loop over sections
    } // End loop over realizations

   // Averaging
    for(int i = 0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            walk_R2_sections[i][j] /= Nrealizations;
        } // End loop over steps (len_sections)
    } // End loop over sections

    // rms
    for(int i=0; i<Nsections; i++){
        for(int j=0; j<len_sections; j++){
            for(int k=0; k<Nrealizations; k++){
                walk_R2_sections_rms[i][j] = (walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j])*(walk_R2_sections_store[k][i][j]-walk_R2_sections[i][j]);
            } // End loop over realizations
            walk_R2_sections_rms[i][j] = sqrt(walk_R2_sections_rms[i][j]/(Nrealizations-1));
        } // End loop over steps (len_sections)
    } // End loop over sections


    //----- Writing to file -----//
    //string outfilenamePrefix;


    ofstream outFile;
    char *filename = new char[100000]; // Possibly a long file name
    sprintf(filename, "Nsteps%i_Nreal%i_pot_exp%f_sigma%f_factor%f_R2_basic_mc", Nrun, Nrealizations, power, sigma, soften);
    outFile.open(filename);
    delete filename;


    for(int i=0; i<Nrun; i++){
        outFile << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl;
        //outFile << std::setprecision(std::numeric_limits<double>::digits10+1) << i+1 << " " << walk_R2[i] << " " << walk_R2_rms[i] << endl; // Is this long printing really neccessary? // Possibly (probably?) gonna use this later.
    }
    outFile.close();

    ofstream poutFile;
    char *pfilename = new char[100000]; // Possibly a long file name
    sprintf(pfilename, "Nsteps%i_Nreal%i_Npart%i_pot_exp%f_sigma%f_factor%f_R2_basic_mc", Nrun, Nrealizations, Nsections, power, sigma, soften);
    poutFile.open(pfilename);
    delete pfilename;

    for(int i=0; i<Nsections; i++){
        poutFile << "--- Section " << i <<" ---" << endl;
        for(int j=0; j<len_sections; j++){
            poutFile << j << " " << walk_R2_sections[i][j] << " " << walk_R2_sections_rms[i][j] << endl;
        } // End loop over steps (len_sections)
    } // End loop over sections
    poutFile.close();

    ofstream coordFile;
    char *cfilename = new char[100000]; // Possibly a long file name
    sprintf(cfilename, "Nsteps%i_Nreal%i_Npart%i_pot_exp%f_sigma%f_factor%f_R2_basic_mc_coords", Nrun, Nrealizations, Nsections, power, sigma, soften);
    coordFile.open(cfilename);
    delete cfilename;

    for(int i=0; i<Nrealizations; i++){
        coordFile << "-- Realization " << i << "--" << endl;
        for(int j=0; j<Nrun; j++){
            coordFile << j << " " << walk_x[i][j] << " " << walk_y[i][j] << endl;
        }
    }
    coordFile.close();

    cout << "acceptance rate: " << acceptance << endl;
    cout << "consecutive_locked: " << consecutive_locked << endl;
    cout << "illegalmovesuggested: " << illegalmovesuggested << endl;
    cout << "e_curr: " << e_curr << endl;
}


vector<int> give_startpoints(int Nrun, int Nsections)
{
    vector<int> startpoints = vector<int>(Nsections);
    int len_sections = Nrun/Nsections;
    for(int i=0; i<Nsections; i++){
        startpoints[i] = i*len_sections;
        cout << "startpoints["<<i<<"]: " << startpoints[i] << endl;
    }
    return startpoints;
}


double give_energy(int thisx, int thisy, int Npots, double sigma, double power, double soften, vector<double> x_centres, vector<double> y_centres){
    // Maybe I should have some sort of cutoff?
    double dx, dy, R2, base;
    double thisenergy;
    double energy = 0;
    //cout << "-------------------" << endl;
    for(int i=0; i<Npots; i++){
        dx      = thisx - x_centres[i];
        dy      = thisy - y_centres[i];
        R2      = dx*dx + dy*dy;
        base    = sigma/R2;
        energy += soften*pow(base, power);
        thisenergy = soften*pow(base, power);
        /*
        cout << "---" << endl;
        cout << "thisx: " << thisx << ", thisy: " << thisy << "; xcentre: " << x_centres[i] << "; ycentre: " << y_centres[i] << endl;
        cout << "distance^2: " << R2 << "; energy contribution: " << thisenergy << "; accumulated energy: " << energy << endl;
        */
    }
    return energy;

}
