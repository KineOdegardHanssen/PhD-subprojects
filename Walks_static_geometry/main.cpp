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

void matrixwalk_hard(int Nblocks, int blocksize, int Nrun, int Nrealizations, int seed);

int main()
{


    int Nrun, Nblocks, Nrealizations, blocksize, seed;//Lx, Ly;
    //Lx = 500;
    //Ly = 500;
    seed = 283;
    Nrun = 20000;
    Nblocks = 3;
    blocksize = 2;
    Nrealizations = 100;

    //srand(seed);

    matrixwalk_hard(Nblocks, blocksize, Nrun, Nrealizations, seed);


    cout << "Hello World!" << endl;
    return 0;
}


void matrixwalk_hard(int Nblocks, int blocksize, int Nrun, int Nrealizations, int seed)
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

    //---- Set up of run ----//
    // Random number distributions
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,Ntotblocks-1);
    std::uniform_int_distribution<int> steps_distr(0,1); // Move: +/-1. Should I be able to keep still too?
    std::uniform_int_distribution<int> dir_distr(0,1);
    //auto location = std::bind ( distribution, generator ); // Couldn't get this to work.

    // Variables
    bool free, moveit;
    //double thisR2;    // Is this a conflict? Should this be an int?
    int starti, startj, indi, indj, nextindi, nextindj, dir, step, thisR2, nx, ny;
    vector<double> walk_R2(Nrun);                                           // Average R^2
    vector<double> walk_R2_rms(Nrun);                                       // RMS R^2
    vector<vector<int>> walk_R2_store(Nrealizations, vector<int>(Nrun)); // For calculation of the standard deviation
    vector<vector<int>> walk_x(Nrealizations, vector<int>(Nrun));        // Useful in case I want to look at the walks for different starting points (post-processing)
    vector<vector<int>> walk_y(Nrealizations, vector<int>(Nrun));

    /*
    vector<int> dir = vector<int>(8); // Direction vector of walks. Inspiration taken from Anders' percwalk.c // Not sure I will keep this. He had flattened the matrix, I think...
    dir[0] = 1;
    dir[1] = 0;
    dir[2] = -1;
    dir[3] = 0;
    dir[4] = 0;
    dir[5] = 1;
    dir[6] = 0;
    dir[7] = -1;
    */
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
            while(!moveit){
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
            //cout << "starti " << starti << "; indi: " << indi << endl;
            //cout << "startj " << startj << "; indj: " << indj << endl;
            //cout << "(indi-starti)*(indi-starti): " << (indi-starti)*(indi-starti) << endl;
            //cout << "(indj-startj)*(indj-startj): " << (indj-startj)*(indj-startj) << endl;
            thisR2 = (indi-starti)*(indi-starti) + (indj-startj)*(indj-startj);
            walk_R2[j] += thisR2;
            walk_R2_store[i][j] = thisR2;
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


}
