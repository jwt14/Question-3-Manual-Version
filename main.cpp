//  HPC_Q3 - Manual version (for reference)
//
//  Created by Jan Witold Tomaszewski CID: 00833865

#include <iostream>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <iomanip>
#include <limits>
#include <fstream>
#include <iterator>
#include "TriMatrix.h"
#include "cblas.h"
using namespace std;

/*----------------------------- Generating Identity and Tri-diagonal Spatial Matrices using help functions ------------------------------------------------------------*/
TriMatrix MakeIdentityMatrix(int N_x){
     TriMatrix Identity(N_x+2);
     for (int i=2; i <=N_x;i++){
                    Identity(i,i) = 1;
                    Identity(i,i-1) = 0;
                    Identity(i,i+1) = 0;
    }
    Identity(1,1) = 1;
    Identity(N_x+1,N_x+1) = 1;
    Identity(1,2) = 0;
    Identity(N_x+1,N_x) = 0;
    return Identity;
}

TriMatrix MakeSpatialOpMatrix(int N_x){
     TriMatrix Spatial(N_x+2);
     for (int i=2; i <=N_x;i++){
                    Spatial(i,i) = -2;
                    Spatial(i,i-1) = 1;
                    Spatial(i,i+1) = 1;
    }
    Spatial(1,1) = 0;
    Spatial(N_x+1,N_x+1) = 0;
    Spatial(1,2) = 0;
    Spatial(N_x+1,N_x) = 0;

    return Spatial;
}

 /*----------------------------- Printing and input-validating functions for both doubles and integers ----------------------------------------------------------------*/
void print_vector(vector<double> U, char const vector_filename[128]){
    ofstream output_file(vector_filename);
    ostream_iterator<double> output_iterator(output_file, "\n");        //Final solution is saved as a text file
    copy(U.begin(), U.end(), output_iterator);
    int m = U.size();
    cout << "The final solution is:" << endl;                           //Final solution is printed to a terminal
    for (int i=0;i<m;++i){
        cout << U[i] << endl;
    }
}

void validating(double &vIn){
    while(1){
        if(cin.fail() || vIn <=0){
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(),'\n');
            cout<<"You have entered wrong input"<<endl;
            cin>>vIn;
        }
    if(!cin.fail())
        break;
    }
}

void validating(int &vIn){
    while(1){
        if(cin.fail() || vIn <= 0){
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(),'\n');
            cout<<"You have entered wrong input"<<endl;
            cin>>vIn;
        }
    if(!cin.fail())
        break;
    }
}


int main() {
    /*----------------------------- Program description ---------------------------------------------------------------------------------------------------------------*/
    //
    //This is HPC Q3 solution- MANUAL VERSION. It provides a solution to a heat equation problem in a form of a temerature vectors.
    //This is modified version that retains cin user input (to check the results in order to validate the quality of the schemes.
    //
    /*----------------------------- Declaring variables and prompting for input with validation -----------------------------------------------------------------------*/
    double L, T, theta;
    int N_x,N_t;

    cout << "\nLength of a domain = "; cin >> L; validating(L);
    cout << "\nNumber of grid points (20, 30 etc): "; cin >> N_x; validating(N_x);      //only even number can be selected!
    double del_x = L/(double(N_x));
    cout << "\nSpatial step size (del_x) = " << del_x << endl;

    cout << "\nTime of simulation = "; cin >> T; validating(T);
    cout << "\nNumber of time steps: "; cin >> N_t; validating(N_t);
    double del_t = T/(double(N_t));
    cout << "\nTime step size (del_x) = " << del_t << endl;

    double alpha;
    cout << "\nThermal conductivity (alpha) = "; cin >> alpha; validating(alpha);
    double nu = alpha*(del_t/pow(del_x,2));                                             //calculating Courant number
    cout << "\nSpecify theta = "; cin >> theta;                                         //asking for Theta parameter

    /*----------------------------- Generating vectors with initial conditions ----------------------------------------------------------------------------------------*/
    vector<double> u_0, u, u_CN;                            //u_0 stores initial heat distribution;
    for(int j=0; j<N_x+1; j++){                             //u stores heat at next full time step;
         u_0.push_back(j*del_x-pow(j*del_x,2));             //u_CN stores heat at the intermediate step for Crank-Nicolson method.
     }

    /*----------------------------- Generating Identity and Spatial triMatrices using help functions ------------------------------------------------------------------*/
    TriMatrix I(N_x+2);
    TriMatrix l(N_x+2);

    l = MakeSpatialOpMatrix(N_x);                       //Lowercase l as uppercase is reserved for bar length!
    I = MakeIdentityMatrix(N_x);

    /*----------------------------- Selecting solver type based on theta value  ---------------------------------------------------------------------------------------*/
    if (theta==0){
        cout << "Selected Forward Euler Method" << endl;
        if (nu >= 0.5){
        cout << "\nCourant number is " << nu << " and scheme is unstable!" << endl << endl;
        }
        else{
        cout << "\nCourant number is " << nu << " and scheme is stable!" << endl << endl;
        }

        TriMatrix A(N_x+2);
        A = I+l*nu;                                            //Generating A matrix with (v, 1-2v, v)

        for(double k=0;k<N_t;++k){
            u = A * u_0;                                    //Implementing for loop with an overloaded vector-matrix multiplication
            u_0 = u;
        }
        print_vector(u,"FEsolution.dat");
    }

    else if (theta==0.5){
        cout << "Selected Crank-Nicolson Method" << endl;
        cout << "\nCourant number is " << nu << " and scheme is unconditionally stable!" << endl << endl;

        TriMatrix A(N_x+2);
        TriMatrix B(N_x+2);
        A = I+l*theta*nu;                                   //Generating A matrix with (v/2, 1-v, v/2)
        B = I-l*theta*nu;                                   //Generating B matrix with (-v/2, 1+v, -v/2)

        for(double k=0;k<N_t;++k){
            u_CN = A * u_0;                                 //Implementing for loop with overloaded multiplication and inversion operators
            u = B / u_CN;                                   //Using Crank-Nicolson two-step method
            u_0 = u;
        }
        print_vector(u,"CNsolution.dat");
    }

    else if (theta==1){
        cout << "Selected Backward Euler Method" << endl;
        cout << "\nCourant number is " << nu << " and scheme is unconditionally stable!" << endl << endl;

        TriMatrix A(N_x+2);
        A = I-l*nu;                                            //Generating A matrix with (-v, 1+v, -v)

        for(double k=0;k<N_t;++k){
            u = A / u_0;                                    //Implementing for loop with an overloaded matrix inversion operator.
            u_0 = u;
        }
        print_vector(u,"BEsolution.dat");
    }

    else{
        cout << "Invalid value of theta. Try again by selecting either 0, 0.5 or 1." << endl;
    }
}
