// Copyright 2020 Maria
#ifndef INCLUDE_MATRIX_HPP_
#define INCLUDE_MATRIX_HPP_

#include <type_traits>
#include <iostream>
#include <limits>
#include <math.h>

template<class T>
class Matrix{
    static_assert(std::is_arithmetic<T>::value, "Not arithmetic type");
    int rows;
    int columns;
    T** m;

public:
    Matrix() {
      rows = 0;
      cols = 0;
      m = nullptr;
    }
    Matrix(int rows, int columns){
        this->rows = rows;
        this->columns = columns;
        m = new T*[rows];
        for (int i = 0; i < rows; i++){
            m[i] = new T[columns];
            for (int j = 0 ; j < columns; j++){
                m[i][j] = 0;
            }
        }
    }
    T* operator[](int index) const{
        return m[index];
    }
    Matrix & operator =(const Matrix& L){
           columns = L.columns;
           rows = L.rows;
           m = new T*[rows];
           for (int i = 0; i < rows; i++){
               m[i] = new T[columns];
           }
           for (int i = 0; i < rows; i++){
               for (int j=0 ; j < columns; j++){
                   m[i][j] = L.m[i][j];
               }
           }
           return *this;
    }
    int Rows() const{
        return rows;
    }
    int Cols() const{
        return columns;
    }
    Matrix operator+(const Matrix& H) const{
        if (columns != H.columns || rows != H.rows) {
            Matrix<T> error;
            return error;
        }
        Matrix<T> sum(rows, columns);
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < columns; j++){
                sum.m[i][j] = m[i][j]+H.m[i][j];
            }
        }
        return sum;
    }
    Matrix operator-(const Matrix& H) const{
        if (columns != H.columns || rows != H.rows) {
            Matrix<T> error;
            return error;
        }
        Matrix<T> raz(rows, columns);
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < columns; j++){
                raz.m[i][j] = m[i][j]-H.m[i][j];
            }
        }
        return raz;
    }
    Matrix operator*(const Matrix& F){
     if (columns != F.rows){
      Matrix<T> error;
      return error;
     }
     Matrix<T> proiz(rows, F.columns);
      for (int i = 0; i < rows; i++){
       for (int j = 0; j < F.columns; j++){
        for (int k = 0; k < columns; k++)
         proiz[i][j] += m[i][k] * F.m[k][j];
       }
      }
      for (int i = 0; i < rows; i++){
       for (int j = 0; j < F.columns; j++){
        if (abs(proiz[i][j])<10*std::numeric_limits<double>::epsilon()){
         proiz[i][j] = 0;
        }
        int w = round(proiz[i][j]);
        if (abs(abs(proiz[i][j])-w)<10*std::numeric_limits<double>::epsilon()){
         proiz[i][j] = round(proiz[i][j]);
        }
       }
      }
     return proiz;
    }
    bool operator != (const Matrix& N) const{
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < columns; j++){
                if (m[i][j] != N[i][j]) {}
                else
                    return false;
            }
        }
        return true;
    }
    bool operator == (const Matrix& O) const{
           for (int i = 0; i < rows; i++){
               for (int j = 0; j < columns; j++){
                   if (m[i][j] == O[i][j]) {}
                   else
                       return false;
               }
           }
           return true;
       }
    ~Matrix(){
        for (int i = 0; i < rows; i++){
            delete[] m[i];
        }
        delete[] m;
    }
    Matrix Inverse() const{
        if (rows != columns){
            Matrix<T> error;
            return error;
        }
        Matrix<T> Turn(rows, 2*columns);
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < columns; j++){
                Turn[i][j] = m[i][j];
            }
        }
        for (int i = 0; i < rows; i++){
            for (int j = columns; j < 2*columns; j++){
                if (i == j-columns) Turn[i][j] = 1;
            }
        }
        double l = 0;
        double n = 0;
        for (int j = 0; j < 2*columns; j++){
            for (int i = 0; i < rows; i++){
                if (i == j){
                    if (Turn[i][j] == 0){
                        int r = i;
                        while (Turn[r][j] == 0){
                            r++;
                        }
                        for (int q = 0; q < 2*columns; q++){
                            double k = Turn[i][q];
                            Turn[i][q] = Turn[r][q];
                            Turn[r][q] = k;
                        }
                    }
                        l = Turn[i][j];
                        for (int k = 0; k < 2*columns; k++){
                            if (Turn[i][k] != 0){
                                Turn[i][k] = Turn[i][k]/l;
                            }
                        }
                }
                if (i > j){
                    n = Turn[i][j]/Turn[j][j];
                    for (int k = 0; k < 2*columns; k++){
                        Turn[i][k] -= Turn[j][k]*n;
                    }
                }
            }
        }
        for (int j = columns-1; j >= 0; j--){
            for (int i = rows-1; i >= 0; i--){
                if (i < j){
                    if (Turn[j][j] != 0){n = Turn[i][j]/Turn[j][j];
                        for (int k = 0; k < 2*columns; k++){
                            Turn[i][k] -= Turn[j][k]*n;
                        }
                    }
                }
            }
        }
        Matrix<T> F(rows, columns);
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < columns; j++){
                F[i][j] = Turn[i][j+columns];
            }
        }
    return F;
    }
    template<class R>
    friend bool operator==(const Matrix<R>& P, const Matrix<R>& Q);
};

template<class T>
bool operator==(const Matrix<double> &P, const Matrix<double> &Q){
  if (std::is_floating_point<T>::value){
    for (int i = 0; i < P.Rows(); i++){
      for (int j = 0; j < P.Cols(); j++){
        if (abs(P[i][j]-Q[i][j])>std::numeric_limits<double>::epsilon()){
          return false;
        }
      }
    }
    return true;
  }
}

#endif // INCLUDE_MATRIX_HPP_

