//  HPC_Q3
//
//  Created by Jan Witold Tomaszewski  CID: 00833865.

#include <iostream>
#include <iomanip>
#include <cstring>
#include <vector>
#include "TriMatrix.h"

using namespace std;

TriMatrix::TriMatrix(const unsigned int pSize)
        : mSize(pSize),
          mZero(0.0) {
    if (pSize == 0) {
        cout << "Matrix size must be > 0" << endl;
        throw;
    }
    mDiag  = new double[mSize];
    mLower = new double[mSize - 1];
    mUpper = new double[mSize - 1];
}

TriMatrix::TriMatrix(const TriMatrix& pSrc)
        : mSize(pSrc.mSize),
          mZero(0.0) {
    mDiag  = new double[mSize];
    mLower = new double[mSize - 1];
    mUpper = new double[mSize - 1];
    memcpy(mDiag,  pSrc.mDiag,  sizeof(double) * mSize);
    memcpy(mLower, pSrc.mLower, sizeof(double) * (mSize - 1));
    memcpy(mUpper, pSrc.mUpper, sizeof(double) * (mSize - 1));
}

TriMatrix::~TriMatrix() {
    delete[] mDiag;
    delete[] mLower;
    delete[] mUpper;
}


TriMatrix& TriMatrix::operator=(const TriMatrix& pSrc) {
    if (mSize != pSrc.mSize) {
        mSize = pSrc.mSize;
        delete[] mDiag;
        delete[] mLower;
        delete[] mUpper;
        mDiag  = new double[mSize];
        mLower = new double[mSize - 1];
        mUpper = new double[mSize - 1];
    }

    memcpy(mDiag,  pSrc.mDiag,  sizeof(double) * mSize);
    memcpy(mLower, pSrc.mLower, sizeof(double) * (mSize - 1));
    memcpy(mUpper, pSrc.mUpper, sizeof(double) * (mSize - 1));
    return *this;
}

double& TriMatrix::operator() (unsigned int pRow, unsigned int pCol) {          //Indexing operator
    if (pRow >= mSize || pCol >= mSize) {                                       //Allows to allocate numbers to a TriMatrix
        cout << "Index out of range." << endl;                                  //Starts with 1 (not with 0) as in usual Mathematic notation
        throw;
    }
    else if (pRow == pCol) {
        return mDiag[pRow-1];
    }
    else if (pCol == pRow-1) {
        return mLower[pRow-1];
    }
    else if (pCol == pRow+1) {
        return mUpper[pRow-1];
    }
    else {
        mZero = 0;
        return mZero;
    }
}

TriMatrix TriMatrix::operator+(const TriMatrix& pSrc) {
    if (mSize != pSrc.mSize) {
        cout << "Matrices are of different sizes." << endl;
        throw;
    }
    TriMatrix result(*this);
    result += pSrc;
    return result;
}

TriMatrix TriMatrix::operator-(const TriMatrix& pSrc) {
    if (mSize != pSrc.mSize) {
        cout << "Matrices are of different sizes." << endl;
        throw;
    }
    TriMatrix result(*this);
    result -= pSrc;
    return result;
}

TriMatrix TriMatrix::operator*(const double& pVal) {
    TriMatrix result(*this);
    result *= pVal;
    return result;
}

TriMatrix TriMatrix::operator/(const double& pVal) {
    TriMatrix result(*this);
    result *= 1.0/pVal;
    return result;
}

TriMatrix& TriMatrix::operator+=(const TriMatrix& pSrc) {
    if (mSize != pSrc.mSize) {
        cout << "Matrices are of different sizes." << endl;
        throw;
    }
    for (unsigned int i = 0; i < mSize; ++i) {
        mDiag[i] += pSrc.mDiag[i];
    }
    for (unsigned int i = 0; i < mSize-1; ++i) {
        mLower[i] += pSrc.mLower[i];
    }
    for (unsigned int i = 0; i < mSize-1; ++i) {
        mUpper[i] += pSrc.mUpper[i];
    }
    return *this;
}

TriMatrix& TriMatrix::operator-=(const TriMatrix& pSrc) {
    if (mSize != pSrc.mSize) {
        cout << "Matrices are of different sizes." << endl;
        throw;
    }
    for (unsigned int i = 0; i < mSize; ++i) {
        mDiag[i] -= pSrc.mDiag[i];
    }
    for (unsigned int i = 0; i < mSize-1; ++i) {
        mLower[i] -= pSrc.mLower[i];
    }
    for (unsigned int i = 0; i < mSize-1; ++i) {
        mUpper[i] -= pSrc.mUpper[i];
    }
    return *this;
}

TriMatrix& TriMatrix::operator*=(const double& pVal) {
    for (unsigned int i = 0; i < mSize; ++i) {
        mDiag[i] *= pVal;
    }
    for (unsigned int i = 0; i < mSize-1; ++i) {
        mLower[i] *= pVal;
    }
    for (unsigned int i = 0; i < mSize-1; ++i) {
        mUpper[i] *= pVal;
    }
    return *this;
}

TriMatrix& TriMatrix::operator/=(const double& pVal) {
    double vTmp = 1.0/pVal;
    for (unsigned int i = 0; i < mSize; ++i) {
        mDiag[i] *= vTmp;
    }
    for (unsigned int i = 0; i < mSize-1; ++i) {
        mLower[i] *= vTmp;
    }
    for (unsigned int i = 0; i < mSize-1; ++i) {
        mUpper[i] *= vTmp;
    }
    return *this;
}

vector<double> TriMatrix::operator* ( vector<double> U)  {              //Vector x Matrix multiplication
    int m = U.size();
    vector<double> u2(m, 0.0);                                          //Vector storing result

    double sum=0;
    for (int i=1;i<m-1;i++){                                            //Performing multiplication apart from 1st and last rows
        sum=0;                                                          //As the value should remain zero at the boundary conditions
        for(int j=1;j<m-1;j++){
            if(i==j){
                sum=sum+mDiag[i]*U[i];
            }
            else if(j==i-1){
                sum=sum+mLower[j]*U[i-1];
            }
            else if(j==i+1){
                sum=sum+mUpper[j]*U[i+1];
            }
        }
        u2[i]=sum;
    }
    return u2;
}

vector<double> TriMatrix::operator/ (vector<double> U)  {               //Matrix ^-1 x Vector operation
    int n = U.size();                                                   //Using Thomas-Algorithm
    vector<double> u2(n, 0.0);
    u2[0]=0.;                                                           //Vector storing results
    u2[n-1]=0.;
    double mDiagTemp[n];                                                //Initiating temporary diagonal matrix
    mDiagTemp[0]=1;

    for(int i=1; i<n; i++){
        double m = mLower[i-1]/mDiagTemp[i-1];
        mDiagTemp[i] = mDiag[i] - m*mUpper[i-1];
        U[i] -= m*U[i-1];
    }
    for(int j=n-2; j >= 0; j--){
        u2[j] = (U[j]- mUpper[j]*u2[j+1])/mDiagTemp[j];
    }
    return u2;
}

void TriMatrix::print() {
    for (unsigned int i = 1; i < mSize; ++i) {
        for (unsigned int j = 1; j < mSize; ++j) {
            cout << setprecision(3) << setw(5) << this->operator()(i,j);
        }
        cout << endl;
    }
    cout << endl;
}
