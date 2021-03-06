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

    srand(seed);

    matrixwalk_hard(Nblocks, blocksize, Nrealizations, Nrun, seed);


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

    for(int i=0; i<Ntotblocks; i++){
        for(int j=0; j<Ntotblocks; j++){
            cout << "walkmat["<<i<<"]["<<j<< "]: " << walkmat[i][j] << endl;
        }
    }

    /*
    int randno;
    for(int i=0; i<50; i++){
        randno = rand() % 100 +1;
        cout << "randno:" << randno << endl;
    }
    */

    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0,Nmat-1);
    //auto location = std::bind ( distribution, generator ); // Couldn't get this to work.

    //---- Set up of run ----//
    bool free;
    int indi, indj;
    vector<int> dir = vector<int>(8); // Direction vector of walks. Inspiration taken from Anders' percwalk.c // Not sure I will keep this. He had flattened the matrix, I think...
    dir[0] = 1;
    dir[1] = 0;
    dir[2] = -1;
    dir[3] = 0;
    dir[4] = 0;
    dir[5] = 1;
    dir[6] = 0;
    dir[7] = -1;

    for(int i=0; i<Nrealizations; i++){
        // Stuff
        // Start the walk by placing the bead in a random position
        // Check that the position is not occupied
        free = false;
        while(!free){
            indi = distribution(generator); //location();//indi = rand() % Nmat; // Random numbers this way (from <random>) is a bit better than rand() from cstdlib, it seems.
            indj = distribution(generator); //location();//indj = rand() % Nmat;
            //cout << "walkmat["<<indi<<"]["<<indj<< "]: " << walkmat[indi][indj] << endl;
            if(walkmat[indi][indj]==0){free=true;} // Only place the bead on a free site
        }
        //cout << "indi: " << indi << "; indj: " << indj << endl; // Works so far
        for(int j=0; j<Nrun; j++){
            // Move the walker
        }
    }


}
