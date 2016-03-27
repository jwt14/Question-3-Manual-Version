//  HPC_Q3
//
//  Created by Jan Witold Tomaszewski  CID: 00833865.

#ifndef CLASS_TRIMATRIX
#define CLASS_TRIMATRIX
#include <vector>
using namespace std;

class TriMatrix {
    public:
        // Custom constructor
        // No default constructor as size zero matrices make no sense
        TriMatrix(const unsigned int n);

        // Copy constructor is critical for the arithmetic operators below
        TriMatrix(const TriMatrix& pSrc);

        // Essential to deallocate the allocated memory
        ~TriMatrix();

        // Overloading operators between Matrices and scalars / vectors / other matrices
        TriMatrix& operator=  (const TriMatrix& pSrc);
        double&    operator() (unsigned int pRow, unsigned int pCol);
        TriMatrix  operator+  (const TriMatrix& pSrc);
        TriMatrix  operator-  (const TriMatrix& pSrc);
        TriMatrix  operator*  (const double&    pVal);
        vector<double> operator* ( vector<double> U);
        vector<double> operator/ ( vector<double> U);
        TriMatrix  operator/  (const double&    pVal);

        TriMatrix& operator+= (const TriMatrix&    pSrc);
        TriMatrix& operator-= (const TriMatrix&    pSrc);
        TriMatrix& operator*= (const double&    pVal);
        TriMatrix& operator/= (const double&    pVal);


        void print();       //Printing function for TriMatrix

    private:
        unsigned int mSize;
        double* mDiag;
        double* mLower;
        double* mUpper;

        double  mZero;
};

#endif
