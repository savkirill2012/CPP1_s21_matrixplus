#ifndef S21MATRIX_H_
#define S21MATRIX_H_H

#include <cmath>
#include <iostream>

class S21Matrix {
 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other);
  ~S21Matrix();

  S21Matrix& operator*=(const S21Matrix& l);
  S21Matrix operator*(const S21Matrix& r);
  S21Matrix& operator*=(const double r);
  S21Matrix operator*(const double r);
  S21Matrix& operator-=(const S21Matrix& r);
  S21Matrix operator-(const S21Matrix& r);
  S21Matrix& operator+=(const S21Matrix& r);
  S21Matrix operator+(const S21Matrix& r);
  S21Matrix& operator=(const S21Matrix& r);
  S21Matrix& operator=(S21Matrix&& r);
  double& operator()(int i, int j);
  bool operator==(const S21Matrix& r);

  bool EqMatrix(const S21Matrix& other);
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);

  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

  double GetRows();
  double GetCols();

  void SetRows(int rows);
  void SetCols(int cols);

 private:
  int rows_, cols_;
  double** matrix_;

  void Clear() noexcept;
  S21Matrix upTriangleMatrix(const S21Matrix& base);
  bool eqMatrixSize(const S21Matrix& other);
};

#endif