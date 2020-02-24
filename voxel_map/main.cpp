#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <cmath>
#include <string>
//#include <cstdio>
#include<stdio.h>
//#include <ostream> // Don't need this?
// This is buggy. Should I not use qt?

using namespace std;
using std::ofstream; using std::string;

vector<int> int_linspace(int framefirst, int framelast, int Nframes);

int main()
{
    int M, N, gridspacing;
    float wall, charge, T, Kangle, Kbond;
    string filename;

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


    ///*
    ifstream inFile;
    inFile.open(filename);
    if (!inFile)
    {
        cerr << "Could not open file :(" << endl;
        exit(1);
    }
    //*/
    int totframes  = 10000000/10000;
    int framefirst = 0;
    int framelast  = totframes;
    int Nall       = M*N;
    int skiplines  = 9;             // If we hit 'ITEM:', skip this many steps...
    vector<int> indices = vector<int>(totframes);
    //indices = int_linspace(framefirst,framelast,totframes);
    //*
    int framespacing = (framelast-framefirst)/(totframes-1); // Double check
    for (int i; i<totframes-1; i++)
    {
        indices[i] = framefirst + i*framespacing;
        cout << "i=" << i << "; value =" << framefirst + i*framespacing << endl;
    }
    //*/

    //std::ofstream << indices;

    cout << "Hello" << ". " << "It is me" << "\n" << totframes << "\n" << indices[0] << "\n" << indices[1] << endl;
    /*
    cout << "Hello" << endl;//indices[1];
    cout << "It is me" << endl;
    cout << totframes << endl;
    cout << indices[0] << endl;
    //cout << indices[1] << endl;*/

    //char filename2[900];

    //filename =
    //sprintf(filename2, 'chaingrid_quadratic_M%iN%i_gridspacing%i_Langevin_wall%.3f_Kangle%.13f_Kbond%.12f_debye_kappa1_debyecutoff3_charge%i_T%i_theta0is180_twofirst_are_fixed.lammpstrj', M, N, gridspacing, wall, Kangle, Kbond, charge, T);

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
