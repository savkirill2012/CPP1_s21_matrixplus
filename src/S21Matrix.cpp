#include "S21Matrix.h"

// constructors

S21Matrix::S21Matrix()
    : rows_(1), cols_(1), matrix_{new double* [rows_] { nullptr }} {
  try {
    matrix_[0] = new double[cols_]{};
  } catch (std::bad_alloc&) {
    Clear();
    throw;
  }
}

S21Matrix::S21Matrix(int rows, int cols)
    : rows_(rows), cols_(cols), matrix_{new double* [rows_] { nullptr }} {
  try {
    for (int i = 0; i < rows_; ++i) {
      matrix_[i] = new double[cols_]{};
    }
  } catch (std::bad_alloc&) {
    Clear();
    throw;
  }
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_), matrix_{new double* [rows_] {
        nullptr
      }} {
  try {
    for (int i = 0; i < rows_; ++i) {
      matrix_[i] = new double[cols_]{};
    }
  } catch (std::bad_alloc&) {
    Clear();
    throw;
  }

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix&& other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.matrix_ = nullptr;
  other.rows_ = 0;
}

S21Matrix::~S21Matrix() { Clear(); }

// operations

S21Matrix& S21Matrix::operator=(const S21Matrix& r) {
  if (this == &r) return *this;

  S21Matrix temp(r.rows_, r.cols_);
  cols_ = temp.cols_;
  std::swap(rows_, temp.rows_);
  std::swap(matrix_, temp.matrix_);

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] = r.matrix_[i][j];
    }
  }

  return *this;
}

S21Matrix& S21Matrix::operator=(S21Matrix&& r) {
  if (this == &r) return *this;

  cols_ = r.cols_;
  std::swap(rows_, r.rows_);
  std::swap(matrix_, r.matrix_);
  return *this;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& r) {
  SumMatrix(r);
  return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix& r) {
  S21Matrix temp = *this;
  temp += r;
  return temp;
}

double& S21Matrix::operator()(int i, int j) {
  if (i >= rows_ || i < 0 || j >= cols_ || j < 0) {
    throw std::out_of_range("Error: out of range");
  } else {
    return matrix_[i][j];
  }
}

bool S21Matrix::operator==(const S21Matrix& r) { return EqMatrix(r); }

S21Matrix& S21Matrix::operator*=(const S21Matrix& r) {
  MulMatrix(r);
  return *this;
}

S21Matrix S21Matrix::operator*(const S21Matrix& r) {
  S21Matrix temp = *this;
  temp *= r;
  return temp;
}

S21Matrix& S21Matrix::operator*=(const double r) {
  MulNumber(r);
  return *this;
}

S21Matrix S21Matrix::operator*(const double r) {
  S21Matrix temp = *this;
  temp *= r;
  return temp;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& r) {
  SubMatrix(r);
  return *this;
}

S21Matrix S21Matrix::operator-(const S21Matrix& r) {
  S21Matrix temp = *this;
  temp -= r;
  return temp;
}

// functions

bool S21Matrix::EqMatrix(const S21Matrix& other) {
  if (eqMatrixSize(other)) {
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        if (fabs(matrix_[i][j] - other.matrix_[i][j]) > 1e-7) {
          return false;
        }
      }
    }
    return true;
  } else {
    return false;
  }
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (eqMatrixSize(other)) {  // work without this-> ?
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        matrix_[i][j] = matrix_[i][j] + other.matrix_[i][j];
      }
    }
  } else {
    throw std::invalid_argument("Matrix have different size");
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (eqMatrixSize(other)) {
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        matrix_[i][j] = matrix_[i][j] - other.matrix_[i][j];
      }
    }
  } else {
    throw std::invalid_argument("Matrix have different size");
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] = matrix_[i][j] * num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ == other.rows_) {
    S21Matrix temp(rows_, other.cols_);
    std::swap(matrix_, temp.matrix_);
    int loc_column = 0;
    while (loc_column < cols_) {
      for (int i = 0; i < rows_; ++i) {
        for (int j = 0, hig = other.cols_; j < hig; ++j) {
          matrix_[i][j] += temp.matrix_[i][loc_column] * other.matrix_[i][j];
        }
      }
      loc_column++;
    }
    cols_ = other.cols_;
  } else {
    throw std::invalid_argument("Matrix isn't square");
  }
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix temp(cols_, rows_);

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      temp(j, i) = matrix_[i][j];
    }
  }

  return temp;
}

double S21Matrix::Determinant() {
  if (rows_ == cols_) {
    double ret_res = 1;
    S21Matrix temp = upTriangleMatrix(*this);
    for (int i = 0; i < rows_; ++i) {
      ret_res *= temp(i, i);
    }
    return ret_res;
  } else {
    throw std::invalid_argument("Matrix isn't square");
  }
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ == cols_) {
    S21Matrix temp(rows_ - 1, cols_ - 1);
    S21Matrix ret_m(rows_, cols_);
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        int loc_row = 0;
        for (int row = 0; row < rows_; ++row) {
          int loc_col = 0;
          for (int col = 0; col < cols_; ++col) {
            if (row != i && col != j) {
              temp.matrix_[loc_row][loc_col] = matrix_[row][col];
              loc_col++;
            }
          }
          if (row != i) {
            loc_row++;
          }
        }
        ret_m(i, j) = pow(-1, i + j) * temp.Determinant();
      }
    }

    return ret_m;
  } else {
    throw std::invalid_argument("Matrix isn't square");
  }
}

S21Matrix S21Matrix::InverseMatrix() {
  double temp = Determinant();
  if (fabs(0.0 - temp) > 1e-7) {
    S21Matrix temp_m = Transpose();
    temp_m = temp_m.CalcComplements();
    temp_m.MulNumber(1 / temp);

    return temp_m;
  } else {
    throw std::runtime_error("Error: determinant is 0");
  }
}

// private func

S21Matrix S21Matrix::upTriangleMatrix(const S21Matrix& base) {
  S21Matrix temp = base;
  int swap_times = 0;

  for (int loc_colum = 0, last_colum = temp.cols_; loc_colum < last_colum;
       ++loc_colum) {
    // search abs biggest val in col and swap rows with curent pos
    int pos_biggest = loc_colum;
    for (int i = loc_colum; i < last_colum; ++i) {
      if (fabs(temp(i, loc_colum)) > fabs(temp(pos_biggest, loc_colum))) {
        pos_biggest = i;
      }
    }

    if (pos_biggest != loc_colum) {
      std::swap(temp.matrix_[loc_colum], temp.matrix_[pos_biggest]);
      swap_times++;
    }

    // if this val == 0 search nearest colm with wall and swap
    if (fabs(temp(loc_colum, loc_colum) - 0) < 1e-7) {
      for (int i = loc_colum; i < last_colum; ++i) {
        if (fabs(temp(loc_colum, i)) > 0) {
          for (int j = 0; i < last_colum; ++j) {
            std::swap(temp(j, loc_colum), temp(j, i));
          }
        }
      }
      swap_times++;
    }

    if (fabs(temp(loc_colum, loc_colum) - 0) < 1e-7) {
      break;
    }

    // to zero all wals in a row
    for (int i = loc_colum + 1; i < last_colum; ++i) {
      double tmp = temp(i, loc_colum);
      for (int j = loc_colum; j < last_colum; ++j) {
        temp(i, j) = temp(i, j) + (-1 * temp(loc_colum, j) * tmp /
                                   temp(loc_colum, loc_colum));
      }
    }
  }

  if (swap_times % 2 == 1) {
    temp(0, 0) *= -1;
  }

  return temp;
}

void S21Matrix::Clear() noexcept {
  for (int i = 0; i < rows_; ++i) {
    delete[] matrix_[i];
  }

  delete[] matrix_;
}

bool S21Matrix::eqMatrixSize(const S21Matrix& other) {
  if (rows_ == other.rows_ && cols_ == other.cols_) {
    return true;
  } else {
    return false;
  }
}

// Gets

double S21Matrix::GetRows() { return rows_; }

double S21Matrix::GetCols() { return cols_; }

void S21Matrix::SetRows(int rows) {
  S21Matrix tmp(rows, cols_);
  std::swap(matrix_, tmp.matrix_);
  int min_row = std::min(rows, rows_);
  std::swap(tmp.rows_, rows_);

  for (int i = 0; i < min_row; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] = tmp.matrix_[i][j];
    }
  }
}

void S21Matrix::SetCols(int cols) {
  S21Matrix tmp(rows_, cols);
  std::swap(matrix_, tmp.matrix_);
  int min_col = std::min(cols, cols_);
  std::swap(tmp.cols_, cols_);

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < min_col; ++j) {
      matrix_[i][j] = tmp.matrix_[i][j];
    }
  }
}