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

void matrixwalk_hard(int Nblocks, int blocksize, int Nrun, int Nrealizations, int Nsections, int seed);
vector<int> give_startpoints(int Nrun, int Nsections);

int main()
{


    int Nrun, Nblocks, Nrealizations, Nsections, blocksize, seed;//Lx, Ly;
    //Lx = 500;
    //Ly = 500;
    seed = 283;
    Nrun = 20000;
    Nblocks = 3;
    blocksize = 2;
    Nsections = 5;
    Nrealizations = 100;

    //srand(seed);

    //vector<int> startpoints = give_startpoints(Nrun, Nsections); // Call this inside function instead
    matrixwalk_hard(Nblocks, blocksize, Nrun, Nrealizations, Nsections, seed);


    cout << "Hello World!" << endl;
    return 0;
}


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
