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
void matrixwalk_hard_easyinput(int sigma, int d, int Nblocks, int Nrun, int Nrealizations, int Nsections, int seed, bool printmatrix);
void matrixmc_hard_easyinput(int sigma, int d, int Nblocks, int Nrun, int Nrealizations, int Nsections, int seed, bool printmatrix);
void matrixmc_potential_easyinput(int intsigma, int intd, int Nblocks, int Nrun, int Nrealizations, int Nsections, int seed, double beta, double sigma, double power, double soften); // Power is expensive


// I could consider removing these functions:
void matrixwalk_hard(int Nblocks, int blocksize, int Nrun, int Nrealizations, int Nsections, int seed);
void matrixmc_hard(int Nblocks, int blocksize, int Nrun, int Nrealizations, int Nsections, int seed);
void matrixmc_potential(int Nblocks, int blocksize, int Nrun, int Nrealizations, int Nsections, int seed, double beta, double sigma, double power, double soften); // Power is expensive

// Auxilliary functions:
vector<int> give_startpoints(int Nrun, int Nsections);
double give_energy(int thisx, int thisy, int Npots, double sigma, double power, double soften, vector<double> x_centres, vector<double> y_centres);

int main()
{
    int intd, intsigma, Nrun, Nblocks, Nrealizations, Nsections, blocksize, seed;//Lx, Ly;
    double beta, sigma, power, soften;
    bool printmatrix;

    // All
    Nrun = 20000;
    Nblocks = 3;
    Nrealizations = 100;

    // Blocky implementation (spacing and obstacles are the same length)
    blocksize = 2;
    Nsections = 5;

    // Flexible hard-block implementation (obstacles and spacing can be of different lengths)
    intd = 10;
    intsigma = 2;
    printmatrix = true;

    // Potential
    beta = 3.0;
    sigma = 2.0;
    power = 6;
    soften = 1;

    // Deprecated. I should remove this
    seed = 283;
    //srand(seed);

    // Function calls:
    //vector<int> startpoints = give_startpoints(Nrun, Nsections); // Call this inside function instead
    //matrixwalk_hard(Nblocks, blocksize, Nrun, Nrealizations, Nsections, seed);
    //matrixmc_hard(Nblocks, blocksize, Nrun, Nrealizations, Nsections, seed);
    //matrixmc_potential(Nblocks, blocksize, Nrun, Nrealizations, Nsections, seed, beta, sigma, power, soften);
    //matrixwalk_hard_easyinput(intsigma, intd, Nblocks, Nrun, Nrealizations, Nsections, seed, printmatrix);
    matrixmc_potential_easyinput(intsigma, intd, Nblocks, Nrun, Nrealizations, Nsections, seed, beta, sigma, power, soften); // Is there any
    //matrixmc_hard_easyinput(intsigma, intd, Nblocks, Nrun, Nrealizations, Nsections, seed, printmatrix);

    cout << "Hello World!" << endl;
    return 0;
}

//---------------- Flexible functions ----------------//
void matrixwalk_hard_easyinput(int sigma, int d, int Nblocks, int Nrun, int Nrealizations, int Nsections, int seed, bool printmatrix)
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
void matrixmc_hard_easyinput(int sigma, int d, int Nblocks, int Nrun, int Nrealizations, int Nsections, int seed, bool printmatrix)
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


void matrixmc_potential_easyinput(int intsigma, int intd, int Nblocks, int Nrun, int Nrealizations, int Nsections, int seed, double beta, double sigma, double power, double soften)
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
                if(nextindi>Nmat || nextindi<-1){
                    nextindi = indi;
                    nextindj = indj;
                    illegalmovesuggested += 1;
                }
            }
            else{
                //cout << "dir y chosen" << endl;
                nextindi = indi;
                nextindj = indj + step;
                if(nextindj>Nmat || nextindj<-1){
                    nextindi = indi;
                    nextindj = indj;
                    illegalmovesuggested += 1;
                }
            }
            e_next    = give_energy(indi, indj, Npots, sigma, power, soften, centres_x, centres_y);
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



//---------------- Not so flexible functions ----------------//
void matrixwalk_hard(int Nblocks, int blocksize, int Nrun, int Nrealizations, int Nsections, int seed)
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

void matrixmc_hard(int Nblocks, int blocksize, int Nrun, int Nrealizations, int Nsections, int seed)
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


void matrixmc_potential(int Nblocks, int blocksize, int Nrun, int Nrealizations, int Nsections, int seed, double beta, double sigma, double power, double soften){
    //---- Matrix set up ----//
    // (to make it the same size as the hard potential one)
    int Ntotblocks, Nmat;
    Nmat = Nblocks*blocksize; // Empty and filled space is of the same size. Maybe a temporary solution?
    Ntotblocks = (2*Nblocks+1)*blocksize;
    vector<vector<int>> walkmat(Ntotblocks, vector<int>(Ntotblocks));

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
                if(nextindi>Ntotblocks || nextindi<-1){
                    nextindi = indi;
                    nextindj = indj;
                    illegalmovesuggested += 1;
                }
            }
            else{
                //cout << "dir y chosen" << endl;
                nextindi = indi;
                nextindj = indj + step;
                if(nextindj>Ntotblocks || nextindj<-1){
                    nextindi = indi;
                    nextindj = indj;
                    illegalmovesuggested += 1;
                }
            }
            e_next    = give_energy(indi, indj, Npots, sigma, power, soften, centres_x, centres_y);
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
